import numpy as np
from cobaya.likelihood import Likelihood

class ClKappaLikelihood(Likelihood):
    # Data and binning matrices (Load from external files if needed)
    bin_edges = np.array([8, 21, 40, 66, 101, 145, 199, 264, 339, 426, 526, 
                           638, 763, 902, 1250, 1300, 1500, 1700, 1800, 
                           2048, 2500, 3000, 5000])
    # Number of alpha parameters is one less than the number of bin edges.
    n_alpha = len(bin_edges) - 1
    def get_requirements(self):
        """
        Tell Cobaya which parameters this likelihood requires.
        This returns a dict with keys "alpha_0", "alpha_1", ..., "alpha_{n_alpha-1}".
        """
        return {f"alpha_{i}": None for i in range(self.__class__.n_alpha)}

    def initialize(self):

        binning_act="/gpfs/fs1/home/r/rbond/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/binning_matrix_act.txt"
        binning_planck="/gpfs/fs1/home/r/rbond/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/binning_matrix_planck.txt"
        binning_spt= "/gpfs/fs1/home/r/rbond/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/binning_spt.txt"
        clkk_act="/gpfs/fs1/home/r/rbond/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/clkk_bandpowers_act.txt"
        clkk_planck="/gpfs/fs1/home/r/rbond/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/clkk_bandpowers_planck.txt"
        clkk_spt="/gpfs/fs1/home/r/rbond/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/clkk_spt.txt"
        clkk_fiducial= "/gpfs/fs1/home/r/rbond/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/clkk_fiducial_2018.txt"
        covariance="/gpfs/fs1/home/r/rbond/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/covmat_actplanckspt3g.txt"
        """Load necessary data and precompute quantities."""
        # Load data and binning matrices (These paths come from the YAML file)
        self.binning_act = np.loadtxt(binning_act)
        self.binning_planck = np.loadtxt(binning_planck)
        self.binning_spt = np.loadtxt(binning_spt)

        self.clkk_act = np.loadtxt(clkk_act)
        self.clkk_planck = np.loadtxt(clkk_planck)
        self.clkk_spt = np.loadtxt(clkk_spt)
        self.clkk_fiducial = np.loadtxt(clkk_fiducial)

        # Load the covariance matrix
        self.covariance = np.loadtxt(covariance)

        # Compute the inverse covariance matrix
        

        # Define the bin edges
        # self.bin_edges = np.array([8, 21, 40, 66, 101, 145, 199, 264, 339, 426, 526, 
        #                            638, 763, 902, 1250, 1300, 1500, 1700, 1800, 
        #                            2048, 2500, 3000, 5000])

        self.bin_edges=np.array([21,40,66,101,145,199,264,339,426,526,638,763,902,1250])
        start = 2
        end = -3
        self.binning_act=self.binning_act[start:end,:]
        nbins_tot_act=self.clkk_act.size
        # Remove trailing bins from ACT part
        sel = np.s_[nbins_tot_act+end:nbins_tot_act]
        self.covariance = np.delete(np.delete(self.covariance,sel,0),sel,1)

        # Remove leading bins from ACT part
        sel = np.s_[:start]
        self.covariance = np.delete(np.delete(self.covariance,sel,0),sel,1)
        self.clkk_act=self.clkk_act[start:end]
        # Compute the full data vector
        self.data = np.concatenate([self.clkk_act, self.clkk_planck, self.clkk_spt])
        nsims=792
        nbins=13
        hartlap_correction = (nsims-nbins-2.)/(nsims-1.)
        self.cov_inv = np.linalg.inv(self.covariance)*hartlap_correction

        # Compute binning matrix
        self.B_matrix = self.top_hat_matrix(np.arange(5000), self.bin_edges)

    def top_hat_matrix(self, L, bin_edges):
        """Creates a top-hat binning matrix."""
        num_bins = len(bin_edges) - 1
        num_L = len(L)
        B = np.zeros((num_bins, num_L))

        for i in range(num_bins):
            L_low, L_high = bin_edges[i], bin_edges[i+1]
            B[i, :] = (L >= L_low) & (L < L_high)

        return B

    def logp(self, **params):
        """Compute the log-likelihood given alpha parameters."""
        # Compute the model
        N_alpha = self.B_matrix.shape[0]
        # Extract the alphas in order: alpha_0, alpha_1, ..., alpha_{N_alpha-1}
        alphas = np.array([params[f"alpha_{i}"] for i in range(N_alpha)])

        model = np.sum(alphas[:, None] * self.B_matrix * self.clkk_fiducial, axis=0)

        # Apply binning matrices
        model_act = self.binning_act @ model[:3000]
        model_planck = self.binning_planck @ model[:3000]
        model_spt = self.binning_spt @ model

        # Concatenate to match the data vector
        model_total = np.concatenate([model_act, model_planck, model_spt])

        residuals = self.data - model_total

        # Compute chi-squared likelihood
        chi2 = residuals.T @ self.cov_inv @ residuals

        return -0.5 * chi2  # Gaussian likelihood