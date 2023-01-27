import numpy as np
import pandas as pd

class DataContainer:
    """
    Allows to load all necessary matrices and vectors.

    :param C: single-cell expression matrix, can be in following format :
                            * a dense 2D numpy.ndarray
    :param a: array of normalized alphas (alpha_tilde). The sum of all alpha_tilde over all cells should be 1 for every patients.
    :param y: array of drug response as survival rate (ratio of cells that survive the given drug treatment).
    :param dC: doubled matrix C. Half of the matrix is associated with label = 1 and other half with label = 0.
    :param dlabs: doubled label array. Half of the array is equal to 1 and other half to 0.
    """
    def __init__(self, S_df, S, W, Y, n_dose=7, dC=None, dlabs=None):
        self.check_input(S, W, Y) # check vector / matrix size
        self.S_df = S_df
        self.S = S
        self.W = W
        self.Y = Y[:,0:n_dose]
        self.Y_pred = Y[:,0:n_dose]
        self.n_dose = n_dose
        self.dS = self.init_double_mat() # doubled matrix
        self.dlabs = self.init_double_labels() # doubled labels
        self.dM7 = self.init_7_mat(self.dS) # 7*2 matrix for 7 drug doses
        self.dL7 = self.init_7_labels(self.dlabs) # 7*2 labels 
        self.Dmat = self.init_7_hotencode(self.dM7, self.S) # get one hot encoder for drug doses
        self.dM7_D7 = np.concatenate([self.dM7, self.Dmat.transpose()]).transpose() # final input matrix

    def check_input(self, S, W, Y):
        check_n_cells = S.shape[1] == W.shape[0] 
        check_n_patients = W.shape[1] == Y.shape[0]
        if check_n_cells & check_n_patients != True:
            raise ValueError('Size of input matrices and vectors do not match')

    def init_7_mat(self, mat):
        return np.tile(mat, (self.n_dose, 1)).transpose()
    
    def init_7_labels(self, lab):
        return np.tile(lab, (1, self.n_dose))[0]

    def init_7_hotencode(self, dM7, mat):
        Dmat = np.zeros((dM7.shape[1], self.n_dose))
        n_cells = mat.shape[1]
        for i in range(Dmat.shape[0]):
            for j in range(Dmat.shape[1]):
                start_1 = j*n_cells*2
                end_1 = start_1 + n_cells*2
                if i >= start_1 :
                    if i < end_1 :
                        Dmat[i,j] = 1
        return Dmat

    def init_double_mat(self):
        return np.concatenate( (self.S, self.S), axis=1).transpose()

    def init_double_labels(self):
        return np.concatenate( ( np.ones(self.S.shape[1]), np.zeros(self.S.shape[1]) ) )
