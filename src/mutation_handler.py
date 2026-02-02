"""
Mutation Handler Module

Handles mutation parsing, validation, and application to sequences.
"""

import re
from typing import List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path

from .sequence_handler import SequenceHandler, VALID_AMINO_ACIDS


@dataclass
class Mutation:
    """Represents a single point mutation."""
    original: str      # Original amino acid (one-letter code)
    position: int      # Position in sequence (1-indexed)
    mutant: str        # Mutant amino acid (one-letter code)
    
    def __post_init__(self):
        """Validate mutation data."""
        if self.original not in VALID_AMINO_ACIDS:
            raise ValueError(f"Invalid original amino acid: {self.original}")
        if self.mutant not in VALID_AMINO_ACIDS:
            raise ValueError(f"Invalid mutant amino acid: {self.mutant}")
        if self.position < 1:
            raise ValueError(f"Position must be >= 1, got {self.position}")
    
    def __str__(self) -> str:
        return f"{self.original}{self.position}{self.mutant}"
    
    def __repr__(self) -> str:
        return f"Mutation({self.original}{self.position}{self.mutant})"
    
    @classmethod
    def from_string(cls, mutation_str: str) -> "Mutation":
        """
        Parse mutation from string format (e.g., 'A123G').
        
        Args:
            mutation_str: Mutation in format [OriginalAA][Position][MutantAA]
            
        Returns:
            Mutation instance
        """
        mutation_str = mutation_str.strip().upper()
        
        # Pattern: letter + number + letter
        pattern = r'^([A-Z])(\d+)([A-Z])$'
        match = re.match(pattern, mutation_str)
        
        if not match:
            raise ValueError(
                f"Invalid mutation format: '{mutation_str}'\n"
                f"Expected format: [OriginalAA][Position][MutantAA]\n"
                f"Example: A123G (Alanine at position 123 to Glycine)"
            )
        
        original, position, mutant = match.groups()
        return cls(original=original, position=int(position), mutant=mutant)


class MutationHandler:
    """Handles mutation operations on sequences."""
    
    def __init__(self, wt_sequence: SequenceHandler):
        """
        Initialize with a wild-type sequence.
        
        Args:
            wt_sequence: Wild-type sequence handler
        """
        self.wt_sequence = wt_sequence
        self.mutations: List[Mutation] = []
    
    def add_mutation(self, mutation: Mutation) -> None:
        """
        Add a mutation after validation.
        
        Args:
            mutation: Mutation to add
        """
        self._validate_mutation(mutation)
        self.mutations.append(mutation)
    
    def add_mutation_from_string(self, mutation_str: str) -> Mutation:
        """
        Parse and add a mutation from string format.
        
        Args:
            mutation_str: Mutation string (e.g., 'A123G')
            
        Returns:
            The parsed Mutation object
        """
        mutation = Mutation.from_string(mutation_str)
        self.add_mutation(mutation)
        return mutation
    
    def _validate_mutation(self, mutation: Mutation) -> None:
        """
        Validate that mutation matches wild-type sequence.
        
        Args:
            mutation: Mutation to validate
        """
        # Check position is within sequence
        if mutation.position > len(self.wt_sequence):
            raise ValueError(
                f"Mutation position {mutation.position} exceeds "
                f"sequence length {len(self.wt_sequence)}"
            )
        
        # Check original residue matches
        wt_residue = self.wt_sequence.get_residue(mutation.position)
        if wt_residue != mutation.original:
            raise ValueError(
                f"Mutation specifies {mutation.original} at position "
                f"{mutation.position}, but wild-type has {wt_residue}"
            )
        
        # Warn if mutation is silent
        if mutation.original == mutation.mutant:
            print(f"Warning: Silent mutation {mutation} (same amino acid)")
    
    def get_mutant_sequence(self, mutation: Optional[Mutation] = None) -> str:
        """
        Generate mutant sequence with applied mutation(s).
        
        Args:
            mutation: Specific mutation to apply. If None, applies all mutations.
            
        Returns:
            Mutant sequence string
        """
        seq_list = list(str(self.wt_sequence))
        
        if mutation:
            mutations_to_apply = [mutation]
        else:
            mutations_to_apply = self.mutations
        
        for mut in mutations_to_apply:
            seq_list[mut.position - 1] = mut.mutant
        
        return "".join(seq_list)
    
    def get_mutant_handler(self, mutation: Mutation) -> SequenceHandler:
        """
        Get a SequenceHandler for the mutant sequence.
        
        Args:
            mutation: Mutation to apply
            
        Returns:
            SequenceHandler for mutant sequence
        """
        mutant_seq = self.get_mutant_sequence(mutation)
        mutant_name = f"{self.wt_sequence.name}_{mutation}"
        return SequenceHandler(mutant_seq, mutant_name)
    
    def clear_mutations(self) -> None:
        """Clear all stored mutations."""
        self.mutations.clear()
    
    def get_mutation_summary(self) -> str:
        """Get a summary of all mutations."""
        if not self.mutations:
            return "No mutations defined"
        
        lines = [f"Wild-type: {self.wt_sequence.name} ({len(self.wt_sequence)} aa)"]
        lines.append(f"Number of mutations: {len(self.mutations)}")
        lines.append("\nMutations:")
        for mut in self.mutations:
            lines.append(f"  - {mut}")
        
        return "\n".join(lines)


def parse_mutation_list(mutation_string: str) -> List[str]:
    """
    Parse a comma or space-separated list of mutations.
    
    Args:
        mutation_string: String containing mutations (e.g., "A123G, V456L")
        
    Returns:
        List of individual mutation strings
    """
    # Split by comma, semicolon, or whitespace
    mutations = re.split(r'[,;\s]+', mutation_string.strip())
    return [m for m in mutations if m]


def create_mutation_directory(
    base_path: Path, 
    wt_name: str, 
    mutation: Mutation
) -> Path:
    """
    Create organized directory structure for a mutation study.
    
    Args:
        base_path: Base output directory
        wt_name: Name of wild-type protein
        mutation: Mutation being studied
        
    Returns:
        Path to mutation-specific directory
    """
    # Create path: base/wt_name/mutation/
    mutation_dir = Path(base_path) / wt_name / str(mutation)
    mutation_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    (mutation_dir / "conformations").mkdir(exist_ok=True)
    (mutation_dir / "aligned").mkdir(exist_ok=True)
    (mutation_dir / "analysis").mkdir(exist_ok=True)
    
    return mutation_dir
