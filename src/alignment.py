"""
Structure Alignment Module

Handles alignment of protein conformations for comparison and analysis.
"""

import numpy as np
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any
from dataclasses import dataclass
import yaml


@dataclass
class AlignmentConfig:
    """Configuration for structure alignment."""
    method: str = "kabsch"           # kabsch, iterative
    selection: str = "backbone"      # backbone, all, ca
    reference: str = "first"         # first, average
    
    @classmethod
    def from_yaml(cls, config_path: Path) -> "AlignmentConfig":
        """Load configuration from YAML file."""
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        align_config = config.get('alignment', {})
        return cls(**{k: v for k, v in align_config.items() if k in cls.__dataclass_fields__})


class StructureAligner:
    """
    Aligns protein conformations using various methods.
    
    Supports:
    - Kabsch algorithm for optimal superposition
    - Iterative alignment with outlier rejection
    - Different atom selections (backbone, CA, all atoms)
    """
    
    def __init__(self, config: Optional[AlignmentConfig] = None):
        """
        Initialize structure aligner.
        
        Args:
            config: Alignment configuration
        """
        self.config = config or AlignmentConfig()
    
    def kabsch_rmsd(
        self,
        P: np.ndarray,
        Q: np.ndarray
    ) -> Tuple[float, np.ndarray, np.ndarray]:
        """
        Calculate RMSD between two structures using Kabsch algorithm.
        
        Args:
            P: First structure coordinates (N x 3)
            Q: Second structure coordinates (N x 3)
            
        Returns:
            Tuple of (RMSD, rotation_matrix, translation_vector)
        """
        # Center structures
        centroid_P = np.mean(P, axis=0)
        centroid_Q = np.mean(Q, axis=0)
        
        P_centered = P - centroid_P
        Q_centered = Q - centroid_Q
        
        # Compute covariance matrix
        H = P_centered.T @ Q_centered
        
        # SVD
        U, S, Vt = np.linalg.svd(H)
        
        # Compute rotation
        R = Vt.T @ U.T
        
        # Handle reflection
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # Apply rotation and compute RMSD
        P_aligned = P_centered @ R
        rmsd = np.sqrt(np.mean(np.sum((P_aligned - Q_centered) ** 2, axis=1)))
        
        # Translation
        t = centroid_Q - centroid_P @ R
        
        return rmsd, R, t
    
    def align_structure(
        self,
        mobile: np.ndarray,
        reference: np.ndarray
    ) -> Tuple[np.ndarray, float]:
        """
        Align mobile structure to reference.
        
        Args:
            mobile: Mobile structure coordinates (N x 3)
            reference: Reference structure coordinates (N x 3)
            
        Returns:
            Tuple of (aligned_coordinates, RMSD)
        """
        rmsd, R, t = self.kabsch_rmsd(mobile, reference)
        aligned = mobile @ R + t
        return aligned, rmsd
    
    def align_trajectory(
        self,
        trajectory: np.ndarray,
        reference_idx: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Align all frames in a trajectory to a reference.
        
        Args:
            trajectory: Trajectory coordinates (n_frames x n_atoms x 3)
            reference_idx: Index of reference frame. If None, uses config setting.
            
        Returns:
            Tuple of (aligned_trajectory, RMSD_array)
        """
        n_frames = trajectory.shape[0]
        
        # Determine reference
        if reference_idx is not None:
            reference = trajectory[reference_idx]
        elif self.config.reference == "first":
            reference = trajectory[0]
        elif self.config.reference == "average":
            reference = np.mean(trajectory, axis=0)
        else:
            reference = trajectory[0]
        
        # Align each frame
        aligned_trajectory = np.zeros_like(trajectory)
        rmsds = np.zeros(n_frames)
        
        for i in range(n_frames):
            aligned_trajectory[i], rmsds[i] = self.align_structure(
                trajectory[i], reference
            )
        
        return aligned_trajectory, rmsds
    
    def compute_average_structure(self, trajectory: np.ndarray) -> np.ndarray:
        """
        Compute average structure from aligned trajectory.
        
        Args:
            trajectory: Aligned trajectory (n_frames x n_atoms x 3)
            
        Returns:
            Average coordinates (n_atoms x 3)
        """
        return np.mean(trajectory, axis=0)
    
    def compute_rmsf(self, trajectory: np.ndarray) -> np.ndarray:
        """
        Compute per-residue RMSF from aligned trajectory.
        
        Args:
            trajectory: Aligned trajectory (n_frames x n_atoms x 3)
            
        Returns:
            RMSF values (n_atoms,)
        """
        average = self.compute_average_structure(trajectory)
        deviations = trajectory - average
        rmsf = np.sqrt(np.mean(np.sum(deviations ** 2, axis=2), axis=0))
        return rmsf


class ConformationManager:
    """
    Manages loading, aligning, and saving conformational ensembles.
    """
    
    def __init__(self, config: Optional[AlignmentConfig] = None):
        """
        Initialize conformation manager.
        
        Args:
            config: Alignment configuration
        """
        self.config = config or AlignmentConfig()
        self.aligner = StructureAligner(config)
        self.conformations: List[np.ndarray] = []
        self.pdb_files: List[Path] = []
    
    def load_conformations_mdtraj(self, pdb_files: List[Path]) -> Any:
        """
        Load conformations using MDTraj.
        
        Args:
            pdb_files: List of PDB file paths
            
        Returns:
            MDTraj trajectory object
        """
        try:
            import mdtraj as md
            
            trajectories = []
            for pdb_file in pdb_files:
                traj = md.load(str(pdb_file))
                trajectories.append(traj)
            
            # Join all trajectories
            combined = md.join(trajectories)
            return combined
            
        except ImportError:
            raise ImportError("MDTraj is required for loading PDB files. Install with: pip install mdtraj")
    
    def load_conformations_biopython(self, pdb_files: List[Path]) -> List[np.ndarray]:
        """
        Load conformations using BioPython.
        
        Args:
            pdb_files: List of PDB file paths
            
        Returns:
            List of coordinate arrays
        """
        try:
            from Bio.PDB import PDBParser
            
            parser = PDBParser(QUIET=True)
            conformations = []
            
            for pdb_file in pdb_files:
                structure = parser.get_structure("protein", str(pdb_file))
                
                # Extract CA coordinates
                coords = []
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if "CA" in residue:
                                coords.append(residue["CA"].get_coord())
                
                conformations.append(np.array(coords))
            
            return conformations
            
        except ImportError:
            raise ImportError("BioPython is required. Install with: pip install biopython")
    
    def align_ensemble(
        self,
        conformations: List[np.ndarray],
        output_dir: Optional[Path] = None
    ) -> Dict[str, Any]:
        """
        Align an ensemble of conformations.
        
        Args:
            conformations: List of coordinate arrays
            output_dir: Directory to save aligned structures
            
        Returns:
            Dictionary with alignment results
        """
        # Stack into trajectory array
        trajectory = np.stack(conformations)
        
        # Align
        aligned_trajectory, rmsds = self.aligner.align_trajectory(trajectory)
        
        # Compute RMSF
        rmsf = self.aligner.compute_rmsf(aligned_trajectory)
        
        # Compute average
        average_structure = self.aligner.compute_average_structure(aligned_trajectory)
        
        results = {
            "n_conformations": len(conformations),
            "n_atoms": conformations[0].shape[0],
            "rmsd_mean": float(np.mean(rmsds)),
            "rmsd_std": float(np.std(rmsds)),
            "rmsd_min": float(np.min(rmsds)),
            "rmsd_max": float(np.max(rmsds)),
            "rmsf_mean": float(np.mean(rmsf)),
            "rmsf_max": float(np.max(rmsf)),
            "aligned_trajectory": aligned_trajectory,
            "rmsds": rmsds,
            "rmsf": rmsf,
            "average_structure": average_structure,
        }
        
        # Save if output directory provided
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            np.save(output_dir / "aligned_trajectory.npy", aligned_trajectory)
            np.save(output_dir / "rmsds.npy", rmsds)
            np.save(output_dir / "rmsf.npy", rmsf)
            np.save(output_dir / "average_structure.npy", average_structure)
        
        return results


def align_conformations_mdtraj(
    pdb_files: List[Path],
    output_dir: Path,
    selection: str = "backbone"
) -> Dict[str, Any]:
    """
    Align conformations using MDTraj (convenience function).
    
    Args:
        pdb_files: List of PDB file paths
        output_dir: Output directory for aligned structures
        selection: Atom selection for alignment
        
    Returns:
        Dictionary with alignment results
    """
    try:
        import mdtraj as md
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load all structures
        print(f"Loading {len(pdb_files)} structures...")
        trajectories = []
        for pdb_file in pdb_files:
            if Path(pdb_file).exists():
                traj = md.load(str(pdb_file))
                trajectories.append(traj)
        
        if not trajectories:
            raise ValueError("No valid PDB files found")
        
        # Combine
        combined = md.join(trajectories)
        
        # Select atoms for alignment
        if selection == "backbone":
            atom_indices = combined.topology.select("backbone")
        elif selection == "ca":
            atom_indices = combined.topology.select("name CA")
        else:
            atom_indices = combined.topology.select("all")
        
        # Align to first frame
        print("Aligning structures...")
        combined.superpose(combined, frame=0, atom_indices=atom_indices)
        
        # Calculate RMSD
        rmsds = md.rmsd(combined, combined, frame=0, atom_indices=atom_indices) * 10  # nm to Angstrom
        
        # Calculate RMSF
        rmsf = md.rmsf(combined, combined, frame=0) * 10  # nm to Angstrom
        
        # Save aligned trajectory
        aligned_traj_path = output_dir / "aligned_ensemble.pdb"
        combined.save(str(aligned_traj_path))
        
        # Save individual frames
        for i in range(combined.n_frames):
            frame_path = output_dir / f"aligned_conf_{i:04d}.pdb"
            combined[i].save(str(frame_path))
        
        # Save metrics
        np.save(output_dir / "rmsds.npy", rmsds)
        np.save(output_dir / "rmsf.npy", rmsf)
        
        results = {
            "n_conformations": combined.n_frames,
            "n_atoms": combined.n_atoms,
            "rmsd_mean": float(np.mean(rmsds)),
            "rmsd_std": float(np.std(rmsds)),
            "aligned_trajectory": str(aligned_traj_path),
            "rmsds": rmsds.tolist(),
            "rmsf": rmsf.tolist(),
        }
        
        print(f"Alignment complete. RMSD: {results['rmsd_mean']:.2f} ± {results['rmsd_std']:.2f} Å")
        
        return results
        
    except ImportError:
        raise ImportError("MDTraj is required. Install with: pip install mdtraj")
