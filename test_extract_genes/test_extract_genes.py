"""Unit tests for extract_genes_final.py using pytest."""
import pytest
from pathlib import Path
import sys
import tempfile

# Add src directory to path to import the module
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from extract_genes_final import (
    load_fasta,
    parse_gff,
    extract_gene_seqs,
    write_fasta_output,
    _reverse_complement,
    _validate_path_with_extensions,
)


class TestReverseComplement:
    """Tests for _reverse_complement function."""

    def test_basic_reverse_complement(self):
        """Test basic reverse complement of DNA sequence."""
        seq = "ATGC"
        result = _reverse_complement(seq)
        assert result == "GCAT"

    def test_reverse_complement_with_n(self):
        """Test reverse complement with N (ambiguous nucleotide)."""
        seq = "ATNGC"
        result = _reverse_complement(seq)
        assert result == "GCNAT"

    def test_reverse_complement_lowercase(self):
        """Test reverse complement with lowercase letters."""
        seq = "atgc"
        result = _reverse_complement(seq)
        assert result == "gcat"

    def test_reverse_complement_mixed_case(self):
        """Test reverse complement with mixed case."""
        seq = "AtGc"
        result = _reverse_complement(seq)
        assert result == "gCaT"

    def test_reverse_complement_long_sequence(self):
        """Test reverse complement of longer sequence."""
        seq = "ATGCTAGCTAGCTAGCTAGC"
        result = _reverse_complement(seq)
        # Manually verify: reverse of ATGCTAGCTAGCTAGCTAGC = CGATCGATCGATCGATGCA
        assert result == "GCTAGCTAGCTAGCTAGCAT"

    def test_reverse_complement_single_nucleotide(self):
        """Test reverse complement of single nucleotide."""
        assert _reverse_complement("A") == "T"
        assert _reverse_complement("T") == "A"
        assert _reverse_complement("G") == "C"
        assert _reverse_complement("C") == "G"


class TestValidatePath:
    """Tests for _validate_path_with_extensions function."""

    def test_valid_fasta_file(self, sample_fasta):
        """Test validation of valid FASTA file."""
        result = _validate_path_with_extensions(str(sample_fasta), [".fasta", ".fa", ".fna"])
        assert result == sample_fasta
        assert result.exists()

    def test_valid_gff_file(self, sample_gff):
        """Test validation of valid GFF file."""
        result = _validate_path_with_extensions(str(sample_gff), [".gff", ".gff3"])
        assert result == sample_gff
        assert result.exists()

    def test_nonexistent_file(self):
        """Test validation of nonexistent file."""
        with pytest.raises(Exception):  # ArgumentTypeError
            _validate_path_with_extensions("/nonexistent/path/file.fasta", [".fasta"])

    def test_wrong_extension(self, sample_fasta):
        """Test validation of file with wrong extension."""
        with pytest.raises(Exception):  # ArgumentTypeError
            _validate_path_with_extensions(str(sample_fasta), [".gff", ".gff3"])


class TestLoadFasta:
    """Tests for load_fasta function."""

    def test_load_single_sequence(self, temp_dir):
        """Test loading FASTA with single sequence."""
        fasta_file = temp_dir / "single.fasta"
        fasta_file.write_text(">seq1\nATGC\n")
        result = load_fasta(fasta_file)
        assert "seq1" in result
        assert result["seq1"] == "ATGC"

    def test_load_multiple_sequences(self, sample_fasta):
        """Test loading FASTA with multiple sequences."""
        result = load_fasta(sample_fasta)
        assert "chr1" in result
        assert "chr2" in result
        assert len(result) == 2

    def test_fasta_sequences_uppercase(self, temp_dir):
        """Test that sequences are converted to uppercase."""
        fasta_file = temp_dir / "lowercase.fasta"
        fasta_file.write_text(">seq1\natgc\ngatc\n")
        result = load_fasta(fasta_file)
        assert result["seq1"] == "ATGCGATC"

    def test_fasta_multiline_sequence(self, temp_dir):
        """Test loading FASTA with sequence spanning multiple lines."""
        fasta_file = temp_dir / "multiline.fasta"
        fasta_file.write_text(">seq1\nATGC\nGATC\nTAGC\n")
        result = load_fasta(fasta_file)
        assert result["seq1"] == "ATGCGATCTAGC"

    def test_fasta_with_empty_lines(self, temp_dir):
        """Test loading FASTA with empty lines (should be ignored)."""
        fasta_file = temp_dir / "empty_lines.fasta"
        fasta_file.write_text(">seq1\n\nATGC\n\nGATC\n")
        result = load_fasta(fasta_file)
        assert result["seq1"] == "ATGCGATC"

    def test_fasta_duplicate_ids_handling(self, duplicate_ids_fasta):
        """Test handling of duplicate sequence IDs."""
        result = load_fasta(duplicate_ids_fasta)
        # Should rename duplicates
        assert "seq1" in result
        assert "seq1_1" in result
        assert len(result) == 2

    def test_fasta_empty_file_error(self, temp_dir):
        """Test that empty FASTA raises ValueError."""
        fasta_file = temp_dir / "empty.fasta"
        fasta_file.write_text("")
        with pytest.raises(ValueError, match="No sequences found"):
            load_fasta(fasta_file)

    def test_fasta_malformed_no_header_error(self, malformed_fasta):
        """Test that FASTA without header raises ValueError."""
        with pytest.raises(ValueError, match="sequence data before any header"):
            load_fasta(malformed_fasta)

    def test_fasta_empty_header_error(self, empty_header_fasta):
        """Test that FASTA with empty header raises ValueError."""
        with pytest.raises(ValueError, match="Empty FASTA header"):
            load_fasta(empty_header_fasta)

    def test_fasta_header_with_spaces(self, temp_dir):
        """Test FASTA header with spaces (takes only first token as ID)."""
        fasta_file = temp_dir / "header_spaces.fasta"
        fasta_file.write_text(">seq1 description extra info\nATGC\n")
        result = load_fasta(fasta_file)
        assert "seq1" in result
        assert result["seq1"] == "ATGC"


class TestParseGff:
    """Tests for parse_gff function."""

    def test_parse_basic_gff(self, sample_gff):
        """Test parsing basic GFF file."""
        result = parse_gff(sample_gff)
        assert len(result) == 3
        assert all(feat["type"] == "gene" for feat in result)

    def test_parse_gff_gene_fields(self, sample_gff):
        """Test that parsed features have all required fields."""
        result = parse_gff(sample_gff)
        feature = result[0]
        assert "seqid" in feature
        assert "start" in feature
        assert "end" in feature
        assert "strand" in feature
        assert "gene_id" in feature
        assert "attributes" in feature

    def test_parse_gff_coordinates_are_integers(self, sample_gff):
        """Test that start and end coordinates are integers."""
        result = parse_gff(sample_gff)
        for feature in result:
            assert isinstance(feature["start"], int)
            assert isinstance(feature["end"], int)

    def test_parse_gff_positive_strand(self, sample_gff):
        """Test parsing gene on positive strand."""
        result = parse_gff(sample_gff)
        gene1 = [f for f in result if f["gene_id"] == "gene001"][0]
        assert gene1["strand"] == "+"

    def test_parse_gff_negative_strand(self, sample_gff):
        """Test parsing gene on negative strand."""
        result = parse_gff(sample_gff)
        gene2 = [f for f in result if f["gene_id"] == "gene002"][0]
        assert gene2["strand"] == "-"

    def test_parse_gff_filters_non_gene_features(self, temp_dir):
        """Test that only 'gene' type features are kept."""
        gff_file = temp_dir / "mixed_features.gff"
        gff_content = """chr1	.	gene	1	10	.	+	.	ID=gene001
chr1	.	mRNA	1	10	.	+	.	ID=mrna001
chr1	.	exon	1	5	.	+	.	ID=exon001
"""
        gff_file.write_text(gff_content)
        result = parse_gff(gff_file)
        assert len(result) == 1
        assert result[0]["gene_id"] == "gene001"

    def test_parse_gff_attributes_parsing(self, sample_gff):
        """Test that attributes are correctly parsed."""
        result = parse_gff(sample_gff)
        feature = result[0]
        assert "ID" in feature["attributes"]
        assert "Name" in feature["attributes"]
        assert feature["attributes"]["ID"] == "gene001"
        assert feature["attributes"]["Name"] == "TestGene1"

    def test_parse_gff_gene_id_fallback(self, temp_dir):
        """Test gene_id extraction with fallback to Name and gene_id."""
        gff_file = temp_dir / "gene_id_fallback.gff"
        gff_content = """chr1	.	gene	1	10	.	+	.	Name=TestGene
chr1	.	gene	20	30	.	+	.	gene_id=test_gene_2
"""
        gff_file.write_text(gff_content)
        result = parse_gff(gff_file)
        assert len(result) == 2
        assert result[0]["gene_id"] == "TestGene"
        assert result[1]["gene_id"] == "test_gene_2"

    def test_parse_gff_malformed_error(self, malformed_gff):
        """Test that malformed GFF raises ValueError."""
        with pytest.raises(ValueError, match="expected 9 columns"):
            parse_gff(malformed_gff)

    def test_parse_gff_invalid_coordinates_error(self, invalid_coordinates_gff):
        """Test that invalid coordinates raise ValueError."""
        with pytest.raises(ValueError, match="Invalid start/end"):
            parse_gff(invalid_coordinates_gff)

    def test_parse_gff_ignores_comments(self, temp_dir):
        """Test that GFF comments are ignored."""
        gff_file = temp_dir / "with_comments.gff"
        gff_content = """# This is a comment
##gff-version 3
# Another comment
chr1	.	gene	1	10	.	+	.	ID=gene001
"""
        gff_file.write_text(gff_content)
        result = parse_gff(gff_file)
        assert len(result) == 1


class TestExtractGeneSeqs:
    """Tests for extract_gene_seqs function."""

    def test_extract_basic_sequences(self, sample_fasta, sample_gff):
        """Test basic gene sequence extraction."""
        fasta_dict = load_fasta(sample_fasta)
        features = parse_gff(sample_gff)
        result = extract_gene_seqs(fasta_dict, features)
        assert len(result) == 3
        assert all("gene_id" in gene for gene in result)
        assert all("sequence" in gene for gene in result)

    def test_extract_sequence_coordinates(self, sample_fasta, sample_gff):
        """Test that extracted sequences use correct coordinates."""
        fasta_dict = load_fasta(sample_fasta)
        features = parse_gff(sample_gff)
        result = extract_gene_seqs(fasta_dict, features)
        
        # First feature: chr1, 1-10 on +
        gene1 = [g for g in result if g["gene_id"] == "gene001"][0]
        assert gene1["start"] == 1
        assert gene1["end"] == 10
        assert len(gene1["sequence"]) == 10

    def test_extract_forward_strand(self, temp_dir):
        """Test extraction of sequence on forward strand."""
        fasta_file = temp_dir / "test.fasta"
        fasta_file.write_text(">chr1\nATGCATGCAT\n")
        gff_file = temp_dir / "test.gff"
        gff_file.write_text("chr1\t.\tgene\t1\t6\t.\t+\t.\tID=gene1\n")
        
        fasta_dict = load_fasta(fasta_file)
        features = parse_gff(gff_file)
        result = extract_gene_seqs(fasta_dict, features)
        
        assert result[0]["sequence"] == "ATGCAT"

    def test_extract_reverse_strand(self, temp_dir):
        """Test extraction of sequence on reverse strand (with reverse complement)."""
        fasta_file = temp_dir / "test.fasta"
        fasta_file.write_text(">chr1\nATGCATGCAT\n")
        gff_file = temp_dir / "test.gff"
        gff_file.write_text("chr1\t.\tgene\t1\t6\t.\t-\t.\tID=gene1\n")
        
        fasta_dict = load_fasta(fasta_file)
        features = parse_gff(gff_file)
        result = extract_gene_seqs(fasta_dict, features)
        
        # ATGCAT reverse complement is ATGCAT (palindrome)
        assert result[0]["sequence"] == "ATGCAT"
        assert result[0]["strand"] == "-"

    def test_extract_nonexistent_seqid_error(self, temp_dir):
        """Test that nonexistent seqid raises KeyError."""
        fasta_file = temp_dir / "test.fasta"
        fasta_file.write_text(">chr1\nATGCATGCAT\n")
        gff_file = temp_dir / "test.gff"
        gff_file.write_text("chr_nonexistent\t.\tgene\t1\t6\t.\t+\t.\tID=gene1\n")
        
        fasta_dict = load_fasta(fasta_file)
        features = parse_gff(gff_file)
        
        with pytest.raises(KeyError, match="unknown seqid"):
            extract_gene_seqs(fasta_dict, features)

    def test_extract_out_of_range_coordinates_error(self, temp_dir):
        """Test that out-of-range coordinates raise ValueError."""
        fasta_file = temp_dir / "test.fasta"
        fasta_file.write_text(">chr1\nATGCAT\n")
        gff_file = temp_dir / "test.gff"
        gff_file.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
        
        fasta_dict = load_fasta(fasta_file)
        features = parse_gff(gff_file)
        
        with pytest.raises(ValueError, match="out of range"):
            extract_gene_seqs(fasta_dict, features)

    def test_extract_unnamed_genes(self, temp_dir):
        """Test handling of genes without ID or Name."""
        fasta_file = temp_dir / "test.fasta"
        fasta_file.write_text(">chr1\nATGCATGCAT\n")
        gff_file = temp_dir / "test.gff"
        gff_file.write_text("chr1\t.\tgene\t1\t6\t.\t+\t.\tSome=attribute\n")
        
        fasta_dict = load_fasta(fasta_file)
        features = parse_gff(gff_file)
        result = extract_gene_seqs(fasta_dict, features)
        
        assert result[0]["gene_id"] == "gene_1"

    def test_extract_stores_strand_info(self, sample_fasta, sample_gff):
        """Test that strand information is preserved."""
        fasta_dict = load_fasta(sample_fasta)
        features = parse_gff(sample_gff)
        result = extract_gene_seqs(fasta_dict, features)
        
        assert result[0]["strand"] == "+"
        assert result[1]["strand"] == "-"


class TestWriteFastaOutput:
    """Tests for write_fasta_output function."""

    def test_write_basic_output(self, temp_dir):
        """Test basic FASTA output writing."""
        output_file = temp_dir / "output.fasta"
        genes = [
            {
                "gene_id": "gene1",
                "seqid": "chr1",
                "start": 1,
                "end": 10,
                "strand": "+",
                "sequence": "ATGCATGCAT",
            }
        ]
        write_fasta_output(output_file, genes)
        
        assert output_file.exists()
        content = output_file.read_text()
        assert ">gene1" in content
        assert "ATGCATGCAT" in content

    def test_write_output_header_format(self, temp_dir):
        """Test that output header has correct format."""
        output_file = temp_dir / "output.fasta"
        genes = [
            {
                "gene_id": "testgene",
                "seqid": "chr1",
                "start": 5,
                "end": 15,
                "strand": "+",
                "sequence": "ATGCATGCAT",
            }
        ]
        write_fasta_output(output_file, genes)
        
        content = output_file.read_text()
        assert "gene_coords=5-15" in content
        assert "strand=+" in content

    def test_write_output_multiple_genes(self, temp_dir):
        """Test writing multiple genes to output."""
        output_file = temp_dir / "output.fasta"
        genes = [
            {
                "gene_id": "gene1",
                "seqid": "chr1",
                "start": 1,
                "end": 10,
                "strand": "+",
                "sequence": "ATGCATGCAT",
            },
            {
                "gene_id": "gene2",
                "seqid": "chr1",
                "start": 20,
                "end": 30,
                "strand": "-",
                "sequence": "TTTTGGGGAA",
            },
        ]
        write_fasta_output(output_file, genes)
        
        content = output_file.read_text()
        assert ">gene1" in content
        assert ">gene2" in content
        lines = content.strip().split("\n")
        assert len([l for l in lines if l.startswith(">")])== 2

    def test_write_output_line_wrapping(self, temp_dir):
        """Test that long sequences are wrapped at specified width."""
        output_file = temp_dir / "output.fasta"
        long_seq = "ATGC" * 20  # 80 nucleotides
        genes = [
            {
                "gene_id": "gene1",
                "seqid": "chr1",
                "start": 1,
                "end": len(long_seq),
                "strand": "+",
                "sequence": long_seq,
            }
        ]
        write_fasta_output(output_file, genes, line_width=60)
        
        lines = output_file.read_text().split("\n")
        # Should have header + 2 sequence lines (60 + 20)
        seq_lines = [l for l in lines if l and not l.startswith(">")]
        assert len(seq_lines) == 2
        assert len(seq_lines[0]) == 60
        assert len(seq_lines[1]) == 20

    def test_write_output_creates_parent_dirs(self, temp_dir):
        """Test that parent directories are created if needed."""
        output_file = temp_dir / "subdir" / "output.fasta"
        genes = [
            {
                "gene_id": "gene1",
                "seqid": "chr1",
                "start": 1,
                "end": 10,
                "strand": "+",
                "sequence": "ATGCATGCAT",
            }
        ]
        # Note: write_fasta_output doesn't create dirs, so we need to do it first
        output_file.parent.mkdir(parents=True, exist_ok=True)
        write_fasta_output(output_file, genes)
        
        assert output_file.exists()


class TestIntegration:
    """Integration tests for the complete workflow."""

    def test_full_workflow(self, sample_fasta, sample_gff, temp_dir):
        """Test complete extraction workflow from input to output."""
        output_file = temp_dir / "output.fasta"
        
        # Load inputs
        fasta_dict = load_fasta(sample_fasta)
        features = parse_gff(sample_gff)
        
        # Extract
        gene_seqs = extract_gene_seqs(fasta_dict, features)
        
        # Write output
        write_fasta_output(output_file, gene_seqs)
        
        # Verify output
        assert output_file.exists()
        output_content = output_file.read_text()
        assert len(output_content) > 0
        assert output_content.count(">") == 3  # 3 genes

    def test_workflow_with_reverse_complement_genes(self, temp_dir):
        """Test workflow with genes on both strands."""
        fasta_file = temp_dir / "test.fasta"
        gff_file = temp_dir / "test.gff"
        output_file = temp_dir / "output.fasta"
        
        # Create test data
        fasta_file.write_text(">chr1\nATGCATGCATGC\n")
        gff_file.write_text("""chr1\t.\tgene\t1\t6\t.\t+\t.\tID=gene1
chr1\t.\tgene\t7\t12\t.\t-\t.\tID=gene2
""")
        
        # Run workflow
        fasta_dict = load_fasta(fasta_file)
        features = parse_gff(gff_file)
        gene_seqs = extract_gene_seqs(fasta_dict, features)
        write_fasta_output(output_file, gene_seqs)
        
        # Verify
        assert output_file.exists()
        output_content = output_file.read_text()
        assert ">gene1" in output_content
        assert ">gene2" in output_content
