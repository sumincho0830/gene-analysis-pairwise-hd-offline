"""Backfill After Preprocess QC for the confirmed-outlier-removal branch."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd

from run_transcriptomics_outlier_removed_diagnostic import (
    write_after_preprocess_qc,
)
from run_transcriptomics_stage1 import PipelineConfig
from run_transcriptomics_variance_filtered import (
    STRATEGIES,
    load_networks,
)
from transcriptomics_qc import build_source_cohort_map


_CONTEXT: dict[str, Any] = {}


def context(run_root: str) -> tuple[PipelineConfig, dict[str, Any], pd.DataFrame, dict[str, str]]:
    if not _CONTEXT:
        config = PipelineConfig(run_root=Path(run_root))
        networks = load_networks(config, STRATEGIES)
        _CONTEXT.update(
            config=config,
            networks={str(network.output_root): network for network in networks},
            clinical=pd.read_csv(config.clinical_path),
            cohort_map=build_source_cohort_map(config.run_root),
        )
    return (
        _CONTEXT["config"],
        _CONTEXT["networks"],
        _CONTEXT["clinical"],
        _CONTEXT["cohort_map"],
    )


def generate_one(run_root: str, output_root: str) -> tuple[str, int]:
    config, networks, clinical, cohort_map = context(run_root)
    network = networks[output_root]
    strategies = [value for value in STRATEGIES if value in network.output_root.parts]
    if len(strategies) != 1:
        raise ValueError(f"Cannot resolve strategy for {network.output_root}")
    paths = write_after_preprocess_qc(
        config,
        network,
        strategies[0],
        clinical,
        cohort_map,
        resume=True,
    )
    return f"{strategies[0]} / {network.network_id}", len(paths)


def pending_network_roots(config: PipelineConfig) -> list[str]:
    pending: list[str] = []
    for network in load_networks(config, STRATEGIES):
        branch = network.output_root / "outlier_removal"
        matrix = branch / "03_preprocess/01_cut_diagnostic/matrix_before_cut.csv.gz"
        complete = branch / "04_after_preprocess/02_data_distribution/pca_scores.csv"
        if matrix.is_file() and not complete.is_file():
            pending.append(str(network.output_root))
    return pending


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    run_root = args.run_root.resolve()
    config = PipelineConfig(run_root=run_root)
    network_roots = pending_network_roots(config)
    print(f"Pending after-preprocess QC networks: {len(network_roots)}", flush=True)
    completed = 0
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {
            executor.submit(generate_one, str(run_root), output_root): output_root
            for output_root in network_roots
        }
        for future in as_completed(futures):
            network_label, path_count = future.result()
            completed += 1
            print(
                f"[{completed}/{len(network_roots)}] {network_label}: {path_count} files",
                flush=True,
            )
    print(f"Completed {completed} networks.", flush=True)


if __name__ == "__main__":
    main()
