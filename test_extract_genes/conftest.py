"""Pytest configuration and shared fixtures for extract_genes tests."""
import pytest
import tempfile
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_fasta(temp_dir):
    """Create a sample FASTA file for testing."""
    fasta_path = temp_dir / "sample.fasta"
    fasta_content = """>chr1
ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>chr2
AAAAAATTTTTTGGGGGGCCCCCC
"""
    fasta_path.write_text(fasta_content)
    return fasta_path


@pytest.fixture
def sample_gff(temp_dir):
    """Create a sample GFF file for testing."""
    gff_path = temp_dir / "sample.gff"
    gff_content = """##gff-version 3
chr1	.	gene	1	10	.	+	.	ID=gene001;Name=TestGene1
chr1	.	gene	20	30	.	-	.	ID=gene002;Name=TestGene2
chr2	.	gene	5	15	.	+	.	ID=gene003;Name=TestGene3
"""
    gff_path.write_text(gff_content)
    return gff_path


@pytest.fixture
def malformed_fasta(temp_dir):
    """Create a malformed FASTA file (no headers)."""
    fasta_path = temp_dir / "malformed.fasta"
    fasta_content = """ATGCTAGCTAGCTAGC
"""
    fasta_path.write_text(fasta_content)
    return fasta_path


@pytest.fixture
def empty_header_fasta(temp_dir):
    """Create a FASTA file with empty header."""
    fasta_path = temp_dir / "empty_header.fasta"
    fasta_content = """>
ATGCTAGCTAGC
"""
    fasta_path.write_text(fasta_content)
    return fasta_path


@pytest.fixture
def duplicate_ids_fasta(temp_dir):
    """Create a FASTA file with duplicate sequence IDs."""
    fasta_path = temp_dir / "duplicate_ids.fasta"
    fasta_content = """>seq1
ATGCTAGCTAGC
>seq1
TTTTGGGGAAAA
"""
    fasta_path.write_text(fasta_content)
    return fasta_path


@pytest.fixture
def malformed_gff(temp_dir):
    """Create a malformed GFF file (wrong column count)."""
    gff_path = temp_dir / "malformed.gff"
    gff_content = """chr1	.	gene	1	10	.	+
"""
    gff_path.write_text(gff_content)
    return gff_path


@pytest.fixture
def invalid_coordinates_gff(temp_dir):
    """Create a GFF file with invalid start/end values."""
    gff_path = temp_dir / "invalid_coords.gff"
    gff_content = """chr1	.	gene	abc	def	.	+	.	ID=gene001
"""
    gff_path.write_text(gff_content)
    return gff_path


@pytest.fixture
def nonexistent_seqid_gff(temp_dir):
    """Create a GFF file with references to non-existent sequence IDs."""
    gff_path = temp_dir / "nonexistent_seqid.gff"
    gff_content = """chr_nonexistent	.	gene	1	10	.	+	.	ID=gene001
"""
    gff_path.write_text(gff_content)
    return gff_path
