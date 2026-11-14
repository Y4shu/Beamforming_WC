# Section VI follow-up: worst-case SOCP SINR vs SNR for K in {30,100,300,exact R} and eps in {2,3,3.1}.
# Run: .venv/bin/python plateau_sweep.py   (prints a table, ~1-2 min)
import numpy as np, cvxpy as cp
N=10; INR=1e3; TRIALS=100
rng=np.random.default_rng(0)
st=lambda t: np.exp(-1j*np.pi*np.arange(N)*np.sin(np.radians(t)))
cg=lambda *s:(rng.standard_normal(s)+1j*rng.standard_normal(s))/np.sqrt(2)
a_bar,a0=st(3.0),st(5.0)
R_in=np.eye(N)+sum(INR*np.outer(st(t),st(t).conj()) for t in [30,50])
def Rhat(snr,K):
    X=cg(N,K)+np.sqrt(10**(snr/10))*np.outer(a0,cg(K))
    for t in [30,50]: X+=np.sqrt(INR)*np.outer(st(t),cg(K))
    return X@X.conj().T/K
def wc(R,eps):
    U=np.linalg.cholesky(R).conj().T; w=cp.Variable(N,complex=True)
    cp.Problem(cp.Minimize(cp.norm(U@w)),[cp.real(a_bar.conj()@w)>=eps*cp.norm(w)+1]).solve()
    return w.value
sinr=lambda w,s:10**(s/10)*abs(w.conj()@a0)**2/np.real(w.conj()@R_in@w)
snrs=np.arange(-10,31,5)
print("SNR:      "+" ".join(f"{s:6d}" for s in snrs))
for K in [30,100,300,None]:
  for eps in [2.0,3.0,3.1]:
    row=[]
    for s in snrs:
        if K is None:  # exact R: no finite-sample error
            R=10**(s/10)*np.outer(a0,a0.conj())+R_in; row.append(10*np.log10(sinr(wc(R,eps),s)))
        else:
            row.append(10*np.log10(np.mean([sinr(wc(Rhat(s,K),eps),s) for _ in range(TRIALS)])))
    print(f"K={str(K or 'inf'):>4} e={eps}: "+" ".join(f"{v:6.1f}" for v in row),flush=True)
