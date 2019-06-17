import pandas as pd
from pathlib import Path
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
from plotnine import *
from collections import Counter
from itertools import combinations
from functools import partial

data = Path("~/Projects/dtpred/data").expanduser()
A = pd.read_msgpack(str(data/"A.msg"))
D = pd.read_msgpack(str(data/"D.msg"))


A_rq = A.groupby(['run','charge'])
A_12 = A_rq.get_group((1,2)).copy()
D_rq = D.groupby(['run','charge'])
D_12 = D_rq.get_group((1,2)).copy()

# how to get rid of the effect of the impact of amino acid counts through the mass?
# project?
# use principal components?
A_12.sequence

(ggplot(A_12) + geom_point(aes("mass", 'dt'), size=1))
LR0 = linregress(A_12.mass, A_12.dt)

def ols_res(df, x,  y):
	M = LinearRegression(copy_X=True, fit_intercept=True)
	x, y = df[[x]], df[[y]]
	M.fit(x,y)
	return y - M.predict(x)

# have to denoise it first, to get some higher quality peptides
A_12['errors'] = ols_res(A_12, 'mass', 'dt')
(ggplot(A_12) + geom_point(aes("mass", 'errors'), size=1))

D_12['errors'] = ols_res(D_12, 'mass', 'dt')
(ggplot(D_12) + geom_point(aes("mass", 'errors'), size=1))


def get_pairs(string, r=2):
	if len(string) > r:
		for i in range(r, len(string)+1):
			yield string[(i-r):i]

def test_get_pairs():
	assert list(get_pairs("ABCAD", 3)) == ['ABC', 'BCA', 'CAD']

def get_counts(df, string_seq_iter=iter, seqstr='sequence'):
	"""Get a DataFrame with counts of string sequences."""
	return pd.DataFrame(Counter(string_seq_iter(s)) for s in df[seqstr]).fillna(0).astype(int)

AAs = get_counts(D_12)
AA_AA = get_counts(D_12, partial(combinations, r=2))
AA_AA_AA = get_counts(D_12, partial(combinations, r=3))
AA_cons = get_counts(D_12, get_pairs)
AAA_cons = get_counts(D_12, partial(get_pairs, r=3))


# shouldn't we model deviations from theoretical mass? Why should these be important???
# because we are about to include the mass???
