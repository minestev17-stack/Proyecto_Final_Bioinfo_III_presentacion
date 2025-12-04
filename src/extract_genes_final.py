"""
extract_genes.py

CLI tool to extract gene sequences from a genome FASTA using a GFF/GFF3 file.

Usage example:
    python extract_genes.py --gff genes.gff --fasta genome.fasta --output genes.fna

"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List


def _validate_path_with_extensions(path_str: str, exts: List[str]) -> Path:
    """Validate that the given path exists and has one of the allowed extensions.

    Raises argparse.ArgumentTypeError on failure.
    """
    p = Path(path_str)
    if not p.exists():
        raise argparse.ArgumentTypeError(f"File not found: {path_str}")
    if p.suffix.lower() not in exts:
        raise argparse.ArgumentTypeError(
            f"File {path_str} must have extension in {exts}"
        )
    return p


def parse_arguments() -> argparse.Namespace:
    """Parse and validate command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes: fasta, gff, output (Path objects).
    """
    parser = argparse.ArgumentParser(
        description="Extract gene sequences from a genome FASTA using a GFF/GFF3 file."
    )

    fasta_type = lambda s: _validate_path_with_extensions(s, [".fa", ".fna", ".fasta"])
    gff_type = lambda s: _validate_path_with_extensions(s, [".gff", ".gff3"])
    # output may not exist yet; validate extension only
    def output_type(s: str) -> Path:
        p = Path(s)
        if p.suffix.lower() not in [".fa", ".fna", ".fasta"]:
            raise argparse.ArgumentTypeError(
                f"Output file must have extension .fa, .fna or .fasta: {s}"
            )
        return p

    parser.add_argument("--fasta", required=True, type=fasta_type, help="Input genome FASTA (.fa/.fna/.fasta)")
    parser.add_argument("--gff", required=True, type=gff_type, help="Input GFF/GFF3 file (.gff/.gff3)")
    parser.add_argument("--output", required=True, type=output_type, help="Output FASTA file (.fa/.fna/.fasta)")

    return parser.parse_args()


def load_fasta(fasta_path: Path) -> Dict[str, str]:
    """Load a FASTA file into a dictionary mapping sequence IDs to sequences.

    Parameters
    ----------
    fasta_path : Path
        Path to FASTA file.

    Returns
    -------
    Dict[str, str]
        Mapping from sequence ID (first token of header) to uppercase sequence string.

    Raises
    ------
    ValueError
        If FASTA is malformed (no headers, empty sequences).
    """
    seqs: Dict[str, List[str]] = {}
    current_id: str | None = None

    with fasta_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                # header
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header in {fasta_path}")
                seq_id = header.split()[0]
                current_id = seq_id
                if current_id in seqs:
                    # append a suffix to make unique
                    suffix = 1
                    new_id = f"{current_id}_{suffix}"
                    while new_id in seqs:
                        suffix += 1
                        new_id = f"{current_id}_{suffix}"
                    current_id = new_id
                seqs[current_id] = []
            else:
                if current_id is None:
                    raise ValueError(f"FASTA file {fasta_path} is malformed: sequence data before any header")
                seqs[current_id].append(line.strip())

    # join and uppercase
    fasta_dict: Dict[str, str] = {k: "".join(v).upper() for k, v in seqs.items()}

    if not fasta_dict:
        raise ValueError(f"No sequences found in FASTA file: {fasta_path}")

    return fasta_dict


def parse_gff(gff_path: Path) -> List[dict]:
    """Parse a GFF/GFF3 file and return a list of gene feature dictionaries.

    Each returned dict contains:
      - seqid (str)
      - source (str)
      - type (str)
      - start (int)
      - end (int)
      - score (str)
      - strand ("+"|"-")
      - phase (str)
      - attributes (dict)

    Only features where type == 'gene' (case-insensitive) are kept.

    Raises
    ------
    ValueError
        If a GFF line is malformed or numeric fields are invalid.
    """
    features: List[dict] = []

    with gff_path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                # allow that some GFFs use spaces, but require 9 standard columns
                raise ValueError(f"Malformed GFF line (expected 9 columns): {line}")
            seqid, source, ftype, start, end, score, strand, phase, attributes_str = parts
            if ftype.lower() != "gene":
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                raise ValueError(f"Invalid start/end in GFF line: {line}")
            # parse attributes (key=value;key2=value2)
            attributes = {}
            for attr in attributes_str.split(";"):
                if not attr:
                    continue
                if "=" in attr:
                    k, v = attr.split("=", 1)
                elif " " in attr:
                    # fallback for some GFFs that use space-separated attributes
                    k, v = attr.split(" ", 1)
                else:
                    k, v = attr, ""
                attributes[k.strip()] = v.strip().strip('"')

            gene_id = attributes.get("ID") or attributes.get("Name") or attributes.get("gene_id")
            features.append(
                {
                    "seqid": seqid,
                    "source": source,
                    "type": ftype,
                    "start": start_i,
                    "end": end_i,
                    "score": score,
                    "strand": strand,
                    "phase": phase,
                    "attributes": attributes,
                    "gene_id": gene_id,
                }
            )

    return features


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Non-ATGC characters are complemented generically (N -> N).
    """
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1]


def extract_gene_seqs(fasta_dict: Dict[str, str], gff_features: List[dict]) -> List[dict]:
    """Extract gene sequences from fasta dictionary using GFF gene features.

    Parameters
    ----------
    fasta_dict : Dict[str,str]
        Mapping seqid -> sequence string.
    gff_features : List[dict]
        List of feature dictionaries returned by `parse_gff`.

    Returns
    -------
    List[dict]
        Each dict contains: gene_id, seqid, start, end, strand, sequence

    Raises
    ------
    KeyError
        If a feature refers to a seqid not present in `fasta_dict`.
    ValueError
        If coordinates are out of range.
    """
    results: List[dict] = []
    unnamed_count = 0

    for feat in gff_features:
        seqid = feat["seqid"]
        if seqid not in fasta_dict:
            raise KeyError(f"GFF feature references unknown seqid '{seqid}'")
        chrom_seq = fasta_dict[seqid]
        start = feat["start"]
        end = feat["end"]
        if start < 1 or end < start or end > len(chrom_seq):
            raise ValueError(
                f"Feature coordinates out of range for {seqid}: {start}-{end} (length {len(chrom_seq)})"
            )
        subseq = chrom_seq[start - 1 : end]
        strand = feat.get("strand", "+")
        if strand == "-":
            subseq = _reverse_complement(subseq)
        gene_id = feat.get("gene_id")
        if not gene_id:
            unnamed_count += 1
            gene_id = f"gene_{unnamed_count}"
        results.append(
            {
                "gene_id": gene_id,
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "sequence": subseq,
            }
        )

    return results


def write_fasta_output(output_path: Path, gene_seqs: List[dict], line_width: int = 60) -> None:
    """Write extracted gene sequences to an output FASTA file.

    Header format:
      >{gene_id}  gene_coords={start}-{end} strand={strand}
    """
    with output_path.open("w", encoding="utf-8") as outfh:
        for gene in gene_seqs:
            header = f">{gene['gene_id']}  gene_coords={gene['start']}-{gene['end']} strand={gene['strand']}"
            outfh.write(header + "\n")
            seq = gene["sequence"]
            for i in range(0, len(seq), line_width):
                outfh.write(seq[i : i + line_width] + "\n")


def main() -> None:
    """Main entrypoint: parse args, run extraction and write output."""
    try:
        args = parse_arguments()
        fasta_dict = load_fasta(args.fasta)
        gff_feats = parse_gff(args.gff)
        gene_seqs = extract_gene_seqs(fasta_dict, gff_feats)
        write_fasta_output(args.output, gene_seqs)
    except Exception as exc:  # top-level catcher to give a clear exit code
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()