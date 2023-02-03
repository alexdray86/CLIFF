.. _detailmethods:

Detailed methods
________________

**Mathematical Modelling of Cell-Type-Specific Drug Sensitivity** 

We assume a cell **c** has two possible outcomes when exposed to a drug **d**: 

.. image:: ../images/cliff_schemes1.png  

If we sum up the probability to survive for each cells we can get the following likelihood function : 

.. math:: \mathcal{L} = \sum_{c=1}^C \delta_{cd} \cdot log P(\delta_{cd} = 1) + ( 1 - \delta_{cd} ) \cdot log ( 1 - P (\delta_{cd}) )

**Problem1**: we cannot directly observe  :math:`\delta_{cd}`. Indeed, drug sensitivity assay is limited to a bulk measurement! 

.. image:: ../images/cliff_schemes2a.png  

**Problem 2**: we have **K** heterogeneous population of cells in the mixture subjected to a drug sensitivity assay, with varying survival rates:  :math:`\pi_{dk} = P(\delta_{cd}^{(k)} = 1)`.

.. image:: ../images/cliff_schemes2b.png  

**CLIFF - Cell fate inference for drug sensitivity deconvolution**

Here we propose to infer the drug sensitivity of individual cell subtypes using only drug sensitivity data observed on mixed cell populations, considering  :math:`K` different cell populations in  :math:`N` mixed samples, with proportions than can be defined as  :math:`\mathbf{w}_n = [w_{n1}, ..., w_{nK}]`, with  :math:`w_{nk} \in [0,1]` and :math:`\sum_{k=1}^K w_{nk} = 1`. Drug sensitivity observed in the bulk samples we wish to demix are represented by a matrix :math:`\mathbf{Y}^{N\times D}`, with :math:`y_{nd} \in [0,1]` being the empirical surviving rate of cells, for sample :math:`n` and drug concentration :math:`d`. Therefore the target matrix we wish to predict is :math:`\hat{\mathbf{Y}}^{N\times D}`, and the predicted values of drug sensitivity is assumed to depend on the underlying drug sensitivity of cell subtypes :math:`\pi_{ndk} \in [0,1]` for patient :math:`n`, drug concentration :math:`d` and cell-subtype :math:`k`, represented as a surviving rate of cells.

The log-likelihood defining our problem, and that CLIFF aims to optimize, is defined by:

.. math:: \mathcal{L} =  \max\limits_{\tilde{\boldsymbol{\beta}}} \sum_{n=1}^N  \sum_{d=1}^D \bigg[ y_{nd} \log\bigg( \sum_{k=1}^K w_{nk} \cdot \pi_{kd} (\tilde{\boldsymbol{\beta}}) \bigg) + ( 1 - y_{nd} ) \log \bigg(1 - \sum_{k=1}^K w_{nk} \cdot \pi_{kd} (\tilde{\boldsymbol{\beta}}) \bigg)  \bigg]

As this log-likelihood expression contains the log of a sum, it is difficult to optimize with respect to :math:`\tilde{\boldsymbol{\beta}}`. Therefore, we use an Expectation-Maximization to iteratively define lower bound functions that we optimize with respect to :math:`\tilde{\boldsymbol{\beta}}`. The EM algorithm rely on iterating over two steps, the E-step and the M-step. During E-step, we define a new lower bound function. During M-step, we find the optimal coefficient for this lower bound function. New coefficient in turns allows to define the next lower bound function, and so on so forth. Here is a summarized description of the E/M steps:

**E-step**

During E-step, we define our lower bound function at step :math:`t`, :math:`\mathcal{L}_B^{(t)}` by defining :math:`q^{(0),(t)}_{ndk}` and :math:`q^{(1),(t)}_{ndk}` using :math:`\pi^{(t-1)}_{dk}` defined at previous step  :

.. math:: q^{(1),(t)}_{ndk}  = \frac{ w_{nk} \cdot \pi^{(t-1)}_{dk} }{\sum_{k=1}^K w_{nk} \cdot \pi^{(t-1)}_{dk} }

.. math:: q^{(0), (t)}_{ndk}  = \frac{ w_{nk} \cdot ( 1 - \pi^{(t-1)}_{dk}) }{\sum_{k=1}^K w_{nk} \cdot  ( 1 - \pi^{(t-1)}_{dk} ) }

This gives us the lower bound function of the log-likelihood function, tight to it whenever :math:`q^{(0),(t)}_{ndk}` and :math:`q^{(1),(t)}_{ndk}` satisfies the two condition above.

**M-step**

During M-step, we optimize the lower bound defined in E-step with respect to :math:`\tilde{\boldsymbol{\beta}}`, considering the variables :math:`\mathbf{q}^{(0), (t)}` and :math:`\mathbf{q}^{(0), (t)}` defined during E-step as constants. Also, the two last term of the lower bound defined above do not depend on :math:`\tilde{\boldsymbol{\beta}}`, thus we can solve the following problem which will produce the same optimal :math:`\tilde{\boldsymbol{\beta}}` coefficients: 

.. math:: \mathcal{L}_B = \max\limits_{\boldsymbol{\tilde{\beta}}} \sum_{d,n}^{D,N} \sum_{k=1}^K q^{(1)}_{ndk} \cdot y_{nd} \cdot  \log( \pi_{kd} (\tilde{\boldsymbol{\beta}}_d ) ) + q^{(0)}_{ndk} \cdot (1 - y_{nd}) \cdot \log( 1 - \pi_{kd} (\tilde{\boldsymbol{\beta}}_d ) )

Interestingly, this equation can be solved with standard libraries using a weighted logistic regression model, where we define weights as being :math:`\boldsymbol{\gamma}^{(1)}_{dk} = \sum_{n=1}^{N} q^{(1)}_{ndk} \cdot y_{nd}` and :math:`\boldsymbol{\gamma}^{(0)}_{dk} = \sum_{n=1}^{N} q^{(0)}_{ndk} \cdot (1 - y_{nd})`.

