import numpy as np
import pandas as pd
import scipy
import random
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from scipy.stats import spearmanr, pearsonr

class Cliff:
    """
    Launch the Expectation-Maximization algorithm by using arrays / matrices provided in the object DC (class DataContainer defined above).

    :param DC: DataContainer object containing input matrix, alphas and target survival rate y
    :param n_iter: integer defining the number of iteration to be used during the EM algorithm 
    """
    def __init__(self, DC, n_iter, param_C=1, n_fold_cv=1, PI_true=None, print_steps=True):
        self.DC = DC # DC is a DataContainer object
        self.n_iter = n_iter # n_iter of EM 
        self.n_fold_cv = n_fold_cv # n fold for cross-validation
        self.param_C = param_C
        self.print_steps = print_steps
        # Variables to be found during the run
        self.Y_pred = None # records betas / y_hat at each iter
        self.Y_true = self.DC.Y
        self.b = None
        self.W = self.DC.W
        self.init_beta0(self.b)
        self.PI = None
        self.PI_true = PI_true
        self.it_loss_neg = 0
        self.loss = []
        # Record of metrics - for bulk prediction and cell-type if provided
        self.bulk_rmse = None ; self.bulk_pcc = None
        self.celltype_rmse = None ; self.celltype_pcc = None
                    
    def run(self):
        b = self.b
        W = self.DC.W
        Y_true = self.Y_true
        N = Y_true.shape[0] ; D = Y_true.shape[1]
        prev_loss = 1e10 ; it_loss_neg = 0
        
        # Iterate over EM steps 
        for x in range(self.n_iter):
            # Expectation-Maximization algorithm
            # E-step to calculate weights used by weighted logistic regression
            weights = self.e_step(self.DC.dM7_D7, W, b, Y_true)
            
            # M-step to update beta coefficient with weighted logistic regression
            b = self.m_step(weights)
            
            # Record results
            self.b = b
            pi = self.get_pi(b)
            self.Y_pred = self.get_y_hats(W, pi)
            
            # Compute loss 
            loss = -(1/(N*D))*np.sum(np.sum(Y_true * np.log(self.Y_pred) + (1 - Y_true) * np.log(1 - self.Y_pred)))
            self.loss.append(loss) # record loss 
            
            # compute metrics between y_hats and y_true 
            self.compute_metrics_bulk()
            
            # if pi_true is provided, compute metrics between pi and pi_true
            if np.any(self.PI_true) != None:  
                self.compute_metrics_celltype()
            
            # Check for convergence
            if self.check_convergence(x, loss, prev_loss):
                break
            prev_loss = loss
        print("Stopping at step " + str(x) + "/" + str(self.n_iter) + ", loss = " + str(loss))
            
    def check_convergence(self, x, loss, prev_loss):
        if prev_loss - loss < 1e-5:
            self.it_loss_neg += 1
            if self.it_loss_neg > 2:
                return(True)
            else:
                return(False)
        elif self.print_steps:
            print("EM step " + str(x) + ", loss = " + str(loss))
            return(False)
    
    def compute_metrics_bulk(self):
        self.bulk_rmse = rmse(self.Y_pred.flatten(), self.Y_true.flatten())
        self.bulk_pcc = pearsonr(self.Y_pred.flatten(), self.Y_true.flatten())[0]
        
    def compute_metrics_celltype(self):
        all_pcc = [] ; all_rmse = []
        for k in range(self.PI.shape[1]):
            all_pcc.append(pearsonr(np.array(self.PI)[:,k], self.PI_true[:,k])[0])
            all_rmse.append(rmse(np.array(self.PI)[:,k], self.PI_true[:,k]))
        self.celltype_pcc = np.mean(all_pcc)
        self.celltype_rmse  = np.mean(all_rmse)
    
    def e_step(self, mat, W, b, Y7):
        weights = np.empty(0)
        for i in range(Y7.shape[1]):
            mat_start = i*2*self.DC.S.shape[1]
            mat_end   = mat_start + self.DC.S.shape[1]
            mat_sub   = mat[mat_start:mat_end,:].T
            pi = np.array([sigmoid(i) for i in np.matmul(b.T, mat_sub)])
            q_1_num = np.array([pi * i for i in W.T]).T
            q_1_den = np.array([pi * i for i in W.T]).T.sum(axis=0)
            q_1 = q_1_num / q_1_den
            q_0_num = np.array([(1 - pi) * i for i in W.T]).T
            q_0_den = np.array([(1 - pi) * i for i in W.T]).T.sum(axis=0)
            q_0 = q_0_num / q_0_den
            y_col = Y7[:,i]
            w_1 = np.array([i * y_col for i in q_1]).sum(axis=1)
            w_0 = np.array([i * (1 - y_col) for i in q_0]).sum(axis=1)
            this_w = np.concatenate((w_1, w_0), axis=None)
            weights = np.concatenate( (weights, this_w), axis=None)
        return weights

    def m_step(self, weights):
        clf = LogisticRegression(max_iter=1e4, solver='liblinear', penalty='l1', C = self.param_C, fit_intercept = False)
        clf = clf.fit(self.DC.dM7_D7, self.DC.dL7.astype('str'), sample_weight = weights)
        return clf.coef_[0]
     
    def get_pi(self, b):
        all_pis = []
        for x in range(self.DC.n_dose): # self.DC.n_dose represents the number of distinct drug doses
            ohe_7 = np.zeros((self.DC.S.shape[1], self.DC.n_dose))
            ohe_7[:,x] = np.ones(self.DC.S.shape[1])
            mat7 = np.concatenate([self.DC.S, ohe_7.T])
            pi = np.array([sigmoid(i) for i in np.matmul(b.T, mat7)])
            all_pis.append(pi)
        self.PI = pd.DataFrame(all_pis, columns = self.DC.S_df.columns)
        return(np.array(self.PI))

    def get_y_hats(self, W, PI):
        Y_hats = np.matmul(W.T, PI.T)
        return(Y_hats)

    def init_beta0(self, b):
        if b == None:
            self.b = np.zeros(self.DC.S.shape[0] + self.DC.n_dose)


class Cliff_CV:
    """
    Launch CLIFF method with or without cross-validation. It also create the container with the required data.

    :param DC: DataContainer object containing input matrix, alphas and target survival rate y
    :param n_iter: integer defining the number of iteration to be used during the EM algorithm 
    """
    def __init__(self, DC, n_iter, param_C=1, PI_true=None, print_steps=True):
        self.DC = DC # DC is a DataContainer object
        self.n_iter = n_iter # n_iter of EM 
        self.param_C = param_C
        self.print_steps = print_steps
        self.PI_true = PI_true
        self.PI_pred = None
        self.cliff_obj = None
        
    def create_cliff_obj(self):
        self.cliff_obj = Cliff(self.DC, n_iter = self.n_iter, param_C = self.param_C, PI_true = self.PI_true, print_steps=self.print_steps)
    
    def run_cliff(self, n_fold_cv):
        self.cliff_obj.run()         
                
    def run_cliff_cv(self, n_fold_cv):
        S = self.cliff_obj.DC.S ; S_df = self.cliff_obj.DC.S_df
        W = self.cliff_obj.DC.W ; Y = self.cliff_obj.Y_true
        PI_true = self.cliff_obj.PI_true
        pi_train = np.zeros((Y.shape[1], S.shape[1]))
        cv_groups = np.random.choice(5, W.shape[1], p=[0.2, 0.2, 0.2, 0.2, 0.2])
        print('Starting Cross-Validation', end=' ')
        for f in range(n_fold_cv):
            sel_test = np.array(cv_groups == f)
            W_train = W[:,~sel_test].copy()
            W_test  = W[:,sel_test].copy()
            Y_train = Y[~sel_test,:].copy()
            Y_test  = Y[sel_test,:].copy()
            
            # Launch EM algorithm and recover predictions
            DC_train = DataContainer(S_df.copy(), S.copy(), W_train, Y_train, n_dose=self.cliff_obj.DC.n_dose)

            # Build EM launcher and define n_iter
            Cliff_obj = Cliff(DC_train, n_iter = self.n_iter, param_C = self.param_C, PI_true = self.PI_true, print_steps=self.print_steps)
            Cliff_obj.run()
            Y_test_pred = Cliff_obj.get_y_hats(W_test, Cliff_obj.PI)
            self.DC.Y_pred[sel_test,:] = Y_test_pred
            pi_train = pi_train + Cliff_obj.PI
            
        self.PI_pred = pi_train / n_fold_cv   
        self.compute_metrics_celltype()
        
    def compute_metrics_celltype(self):
        all_pcc = [] ; all_rmse = []
        for k in range(self.PI_pred.shape[1]):
            all_pcc.append(pearsonr(np.array(self.PI_pred)[:,k], self.PI_true[:,k])[0])
            all_rmse.append(rmse(np.array(self.PI_pred)[:,k], self.PI_true[:,k]))
        self.celltype_pcc = np.mean(all_pcc)
        self.celltype_rmse  = np.mean(all_rmse)
