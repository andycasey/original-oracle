
import numpy as np
from scipy import stats, optimize as op


def line_lnlike(theta, x, y, yerr):

    m, b, Po, Yo, Vo = theta

    line_model = m * x + b
    line_ivar = 1.0/(yerr**2)
    outlier_model = Yo
    outlier_ivar = 1.0/(yerr**2 + Vo)

    model_likelihood = -0.5 * ((y - line_model)**2 * line_ivar - np.log(line_ivar))
    outlier_likelihood = -0.5 * ((y - outlier_model)**2 * outlier_ivar - np.log(outlier_ivar))

    return np.sum(np.logaddexp(
        np.log(1 - Po) + model_likelihood,
        np.log(Po) + outlier_likelihood
    ))


def line_lnprior(theta):

    m, b, Po, Yo, Vo = theta
    if (1 > Po > 0):
        return 0
    return -np.inf
    

def line_lnprob(theta, x, y, yerr):
    return line_lnprior(theta) + line_lnlike(theta, x, y, yerr)
    


def fit(x, y, y_uncertainty=None, outliers=True, full_output=False):

    if y_uncertainty is None or not np.all(np.isfinite(y_uncertainty)):
        if full_output:
            return stats.linregress(x=x, y=y)
        return stats.linregress(x=x, y=y)[0]

    A = np.vstack((np.ones_like(x), x)).T
    C = np.diag(y_uncertainty * y_uncertainty)
    cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
    b, m = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))

    # Should we be modelling outliers?
    if outliers:
        theta_init = m, b, 0.0, np.median(y), 1.
        theta_opt = op.fmin(lambda *args: -line_lnprob(*args), theta_init,
            args=(x, y, y_uncertainty), disp=False)

        if full_output:
            return theta_opt
        return theta_opt[0]

    if full_output:
        return m, b
    return m