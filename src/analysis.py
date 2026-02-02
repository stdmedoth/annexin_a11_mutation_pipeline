"""
Conformational Analysis Module

Tools for analyzing protein conformational ensembles.
"""

import numpy as np
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import json


@dataclass
class AnalysisConfig:
    """Configuration for analysis."""
    rmsd_enabled: bool = True
    rmsf_enabled: bool = True
    rg_enabled: bool = True
    secondary_structure_enabled: bool = True
    contact_map_enabled: bool = True
    pca_enabled: bool = True
    contact_cutoff: float = 8.0
    pca_components: int = 10


class ConformationalAnalyzer:
    """
    Comprehensive analysis of protein conformational ensembles.
    """
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        """
        Initialize analyzer.
        
        Args:
            config: Analysis configuration
        """
        self.config = config or AnalysisConfig()
    
    def compute_radius_of_gyration(self, coords: np.ndarray) -> float:
        """
        Compute radius of gyration for a structure.
        
        Args:
            coords: Atomic coordinates (N x 3)
            
        Returns:
            Radius of gyration in Angstroms
        """
        centroid = np.mean(coords, axis=0)
        distances = np.sqrt(np.sum((coords - centroid) ** 2, axis=1))
        rg = np.sqrt(np.mean(distances ** 2))
        return rg
    
    def compute_rg_trajectory(self, trajectory: np.ndarray) -> np.ndarray:
        """
        Compute radius of gyration for each frame in trajectory.
        
        Args:
            trajectory: Trajectory (n_frames x n_atoms x 3)
            
        Returns:
            Array of Rg values
        """
        rg_values = np.array([
            self.compute_radius_of_gyration(frame) 
            for frame in trajectory
        ])
        return rg_values
    
    def compute_contact_map(
        self,
        coords: np.ndarray,
        cutoff: Optional[float] = None
    ) -> np.ndarray:
        """
        Compute residue contact map.
        
        Args:
            coords: CA coordinates (n_residues x 3)
            cutoff: Distance cutoff for contacts (Angstroms)
            
        Returns:
            Contact map (n_residues x n_residues)
        """
        cutoff = cutoff or self.config.contact_cutoff
        n_residues = coords.shape[0]
        
        # Compute pairwise distances
        distances = np.zeros((n_residues, n_residues))
        for i in range(n_residues):
            for j in range(i, n_residues):
                dist = np.linalg.norm(coords[i] - coords[j])
                distances[i, j] = dist
                distances[j, i] = dist
        
        # Create contact map
        contact_map = (distances < cutoff).astype(float)
        
        return contact_map
    
    def compute_average_contact_map(
        self,
        trajectory: np.ndarray,
        cutoff: Optional[float] = None
    ) -> np.ndarray:
        """
        Compute average contact map over trajectory.
        
        Args:
            trajectory: Trajectory (n_frames x n_residues x 3)
            cutoff: Distance cutoff
            
        Returns:
            Average contact probability map
        """
        contact_maps = np.array([
            self.compute_contact_map(frame, cutoff)
            for frame in trajectory
        ])
        return np.mean(contact_maps, axis=0)
    
    def perform_pca(
        self,
        trajectory: np.ndarray,
        n_components: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        Perform Principal Component Analysis on trajectory.
        
        Args:
            trajectory: Aligned trajectory (n_frames x n_atoms x 3)
            n_components: Number of components to compute
            
        Returns:
            Dictionary with PCA results
        """
        n_components = n_components or self.config.pca_components
        n_frames = trajectory.shape[0]
        
        # Flatten coordinates
        flattened = trajectory.reshape(n_frames, -1)
        
        # Center data
        mean_structure = np.mean(flattened, axis=0)
        centered = flattened - mean_structure
        
        # Compute covariance matrix
        cov_matrix = np.cov(centered.T)
        
        # Eigendecomposition
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
        
        # Sort by eigenvalue (descending)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        
        # Keep top components
        n_components = min(n_components, len(eigenvalues))
        eigenvalues = eigenvalues[:n_components]
        eigenvectors = eigenvectors[:, :n_components]
        
        # Project data
        projections = centered @ eigenvectors
        
        # Explained variance
        total_variance = np.sum(eigenvalues)
        explained_variance = eigenvalues / total_variance
        cumulative_variance = np.cumsum(explained_variance)
        
        return {
            "eigenvalues": eigenvalues,
            "eigenvectors": eigenvectors,
            "projections": projections,
            "explained_variance": explained_variance,
            "cumulative_variance": cumulative_variance,
            "mean_structure": mean_structure.reshape(-1, 3),
            "n_components": n_components,
        }
    
    def compute_free_energy_landscape(
        self,
        projections: np.ndarray,
        pc1: int = 0,
        pc2: int = 1,
        bins: int = 50,
        temperature: float = 300.0
    ) -> Dict[str, Any]:
        """
        Compute 2D free energy landscape from PCA projections.
        
        Args:
            projections: PCA projections (n_frames x n_components)
            pc1: First principal component index
            pc2: Second principal component index
            bins: Number of bins for histogram
            temperature: Temperature in Kelvin
            
        Returns:
            Dictionary with free energy landscape
        """
        kB = 0.001987  # kcal/(mol·K)
        
        # Extract projections for selected PCs
        x = projections[:, pc1]
        y = projections[:, pc2]
        
        # 2D histogram
        hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
        
        # Convert to probability
        prob = hist / np.sum(hist)
        
        # Handle zeros
        prob[prob == 0] = 1e-10
        
        # Convert to free energy
        free_energy = -kB * temperature * np.log(prob)
        free_energy -= np.min(free_energy)  # Set minimum to zero
        
        # Get bin centers
        x_centers = (xedges[:-1] + xedges[1:]) / 2
        y_centers = (yedges[:-1] + yedges[1:]) / 2
        
        return {
            "free_energy": free_energy,
            "x_centers": x_centers,
            "y_centers": y_centers,
            "pc1": pc1,
            "pc2": pc2,
            "temperature": temperature,
            "x_range": (float(x.min()), float(x.max())),
            "y_range": (float(y.min()), float(y.max())),
        }
    
    def full_analysis(
        self,
        trajectory: np.ndarray,
        output_dir: Optional[Path] = None
    ) -> Dict[str, Any]:
        """
        Perform full conformational analysis.
        
        Args:
            trajectory: Aligned trajectory (n_frames x n_atoms x 3)
            output_dir: Directory to save results
            
        Returns:
            Dictionary with all analysis results
        """
        results = {
            "n_frames": trajectory.shape[0],
            "n_atoms": trajectory.shape[1],
        }
        
        # Radius of gyration
        if self.config.rg_enabled:
            rg_values = self.compute_rg_trajectory(trajectory)
            results["rg"] = {
                "values": rg_values.tolist(),
                "mean": float(np.mean(rg_values)),
                "std": float(np.std(rg_values)),
                "min": float(np.min(rg_values)),
                "max": float(np.max(rg_values)),
            }
        
        # Contact map
        if self.config.contact_map_enabled:
            avg_contacts = self.compute_average_contact_map(trajectory)
            results["contact_map"] = {
                "average": avg_contacts.tolist(),
                "n_contacts_mean": float(np.sum(avg_contacts) / 2),
            }
        
        # PCA
        if self.config.pca_enabled:
            pca_results = self.perform_pca(trajectory)
            
            # Free energy landscape
            fel = self.compute_free_energy_landscape(pca_results["projections"])
            
            results["pca"] = {
                "explained_variance": pca_results["explained_variance"].tolist(),
                "cumulative_variance": pca_results["cumulative_variance"].tolist(),
                "projections": pca_results["projections"].tolist(),
            }
            results["free_energy_landscape"] = {
                "x_range": fel["x_range"],
                "y_range": fel["y_range"],
                "energy_min": float(np.min(fel["free_energy"])),
                "energy_max": float(np.max(fel["free_energy"][np.isfinite(fel["free_energy"])])),
            }
        
        # Save results
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Save as JSON (excluding large arrays)
            summary = {k: v for k, v in results.items() 
                      if not isinstance(v, np.ndarray)}
            with open(output_dir / "analysis_summary.json", 'w') as f:
                json.dump(summary, f, indent=2, default=str)
            
            # Save numpy arrays
            if "rg" in results:
                np.save(output_dir / "rg_values.npy", results["rg"]["values"])
            if "pca" in results:
                np.save(output_dir / "pca_projections.npy", pca_results["projections"])
                np.save(output_dir / "pca_eigenvalues.npy", pca_results["eigenvalues"])
        
        return results


class ComparativeAnalyzer:
    """
    Compare conformational ensembles between WT and mutant.
    """
    
    def __init__(self):
        self.analyzer = ConformationalAnalyzer()
    
    def compare_ensembles(
        self,
        wt_trajectory: np.ndarray,
        mutant_trajectory: np.ndarray,
        output_dir: Optional[Path] = None
    ) -> Dict[str, Any]:
        """
        Compare WT and mutant conformational ensembles.
        
        Args:
            wt_trajectory: WT trajectory (n_frames x n_atoms x 3)
            mutant_trajectory: Mutant trajectory (n_frames x n_atoms x 3)
            output_dir: Output directory
            
        Returns:
            Comparison results
        """
        # Analyze both
        wt_analysis = self.analyzer.full_analysis(wt_trajectory)
        mutant_analysis = self.analyzer.full_analysis(mutant_trajectory)
        
        comparison = {
            "wt": wt_analysis,
            "mutant": mutant_analysis,
            "differences": {},
        }
        
        # Compare Rg
        if "rg" in wt_analysis and "rg" in mutant_analysis:
            comparison["differences"]["rg"] = {
                "delta_mean": mutant_analysis["rg"]["mean"] - wt_analysis["rg"]["mean"],
                "wt_mean": wt_analysis["rg"]["mean"],
                "mutant_mean": mutant_analysis["rg"]["mean"],
            }
        
        # Save comparison
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            with open(output_dir / "comparison_summary.json", 'w') as f:
                json.dump(comparison, f, indent=2, default=str)
        
        return comparison


# Visualization utilities
def create_analysis_plots(
    results: Dict[str, Any],
    output_dir: Path
) -> List[Path]:
    """
    Create visualization plots for analysis results.
    
    Args:
        results: Analysis results dictionary
        output_dir: Output directory for plots
        
    Returns:
        List of created plot paths
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print("Matplotlib and seaborn required for plotting")
        return []
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_paths = []
    
    # Plot Rg distribution
    if "rg" in results:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(results["rg"]["values"], bins=30, edgecolor='black', alpha=0.7)
        ax.axvline(results["rg"]["mean"], color='red', linestyle='--', 
                   label=f'Mean: {results["rg"]["mean"]:.2f} Å')
        ax.set_xlabel("Radius of Gyration (Å)")
        ax.set_ylabel("Count")
        ax.set_title("Radius of Gyration Distribution")
        ax.legend()
        path = output_dir / "rg_distribution.png"
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        plot_paths.append(path)
    
    # Plot PCA variance
    if "pca" in results:
        fig, ax = plt.subplots(figsize=(8, 6))
        explained_var = results["pca"]["explained_variance"]
        cumulative_var = results["pca"]["cumulative_variance"]
        
        x = range(1, len(explained_var) + 1)
        ax.bar(x, explained_var, alpha=0.7, label='Individual')
        ax.plot(x, cumulative_var, 'ro-', label='Cumulative')
        ax.set_xlabel("Principal Component")
        ax.set_ylabel("Explained Variance")
        ax.set_title("PCA Variance Explained")
        ax.legend()
        path = output_dir / "pca_variance.png"
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        plot_paths.append(path)
    
    # Plot contact map
    if "contact_map" in results:
        fig, ax = plt.subplots(figsize=(10, 8))
        contact_map = np.array(results["contact_map"]["average"])
        sns.heatmap(contact_map, cmap='YlOrRd', ax=ax)
        ax.set_xlabel("Residue")
        ax.set_ylabel("Residue")
        ax.set_title("Average Contact Map")
        path = output_dir / "contact_map.png"
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        plot_paths.append(path)
    
    return plot_paths
