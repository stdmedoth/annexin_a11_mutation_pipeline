"""
Sequence Handler Module

Handles protein sequence input, validation, and manipulation.
"""

import re
from typing import Optional, Tuple
from pathlib import Path


# Standard amino acid one-letter codes
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


class SequenceHandler:
    """Handles protein sequence operations."""
    
    def __init__(self, sequence: str, name: str = "protein"):
        """
        Initialize with a protein sequence.
        
        Args:
            sequence: Amino acid sequence (one-letter code)
            name: Name identifier for the protein
        """
        self.raw_sequence = sequence
        self.sequence = self._clean_sequence(sequence)
        self.name = name
        self._validate()
    
    def _clean_sequence(self, sequence: str) -> str:
        """Remove whitespace and convert to uppercase."""
        return re.sub(r'\s+', '', sequence.upper())
    
    def _validate(self) -> None:
        """Validate the sequence contains only valid amino acids."""
        invalid_chars = set(self.sequence) - VALID_AMINO_ACIDS
        if invalid_chars:
            raise ValueError(
                f"Invalid amino acid(s) found: {', '.join(invalid_chars)}\n"
                f"Valid amino acids are: {''.join(sorted(VALID_AMINO_ACIDS))}"
            )
        if len(self.sequence) == 0:
            raise ValueError("Sequence cannot be empty")
    
    def __len__(self) -> int:
        return len(self.sequence)
    
    def __str__(self) -> str:
        return self.sequence
    
    def __repr__(self) -> str:
        return f"SequenceHandler(name='{self.name}', length={len(self)})"
    
    def get_residue(self, position: int) -> str:
        """
        Get residue at a specific position (1-indexed).
        
        Args:
            position: Position in sequence (1-indexed)
            
        Returns:
            Single letter amino acid code
        """
        if position < 1 or position > len(self.sequence):
            raise ValueError(
                f"Position {position} out of range. "
                f"Valid range: 1-{len(self.sequence)}"
            )
        return self.sequence[position - 1]
    
    def to_fasta(self, line_length: int = 60) -> str:
        """
        Convert sequence to FASTA format.
        
        Args:
            line_length: Number of characters per line
            
        Returns:
            FASTA formatted string
        """
        lines = [f">{self.name}"]
        for i in range(0, len(self.sequence), line_length):
            lines.append(self.sequence[i:i + line_length])
        return "\n".join(lines)
    
    def save_fasta(self, filepath: Path) -> None:
        """Save sequence to FASTA file."""
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, 'w') as f:
            f.write(self.to_fasta())
    
    @classmethod
    def from_fasta(cls, filepath: Path) -> "SequenceHandler":
        """
        Load sequence from FASTA file.
        
        Args:
            filepath: Path to FASTA file
            
        Returns:
            SequenceHandler instance
        """
        filepath = Path(filepath)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        name = "protein"
        sequence_lines = []
        
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
            elif line:
                sequence_lines.append(line)
        
        sequence = "".join(sequence_lines)
        return cls(sequence, name)


def validate_sequence_input(sequence: str) -> Tuple[bool, Optional[str]]:
    """
    Validate a sequence input string.
    
    Args:
        sequence: Input sequence string
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        SequenceHandler(sequence)
        return True, None
    except ValueError as e:
        return False, str(e)


def get_sequence_stats(sequence: str) -> dict:
    """
    Calculate basic statistics for a sequence.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        Dictionary with sequence statistics
    """
    handler = SequenceHandler(sequence)
    seq = handler.sequence
    
    # Count amino acids
    aa_counts = {aa: seq.count(aa) for aa in VALID_AMINO_ACIDS}
    
    # Calculate properties
    stats = {
        "length": len(seq),
        "amino_acid_counts": aa_counts,
        "molecular_weight_approx": len(seq) * 110,  # Rough average
        "charged_residues": sum(aa_counts.get(aa, 0) for aa in "DEKRH"),
        "hydrophobic_residues": sum(aa_counts.get(aa, 0) for aa in "AVILMFYW"),
        "polar_residues": sum(aa_counts.get(aa, 0) for aa in "STNQ"),
    }
    
    return stats
