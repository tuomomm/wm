from neuron import h
h.load_file("nqs.hoc")
import numpy

NQS = h.NQS
nqsdel = h.nqsdel

# converts 2D numpy array into an NQS
#  NB: each row of the numpy array will be a Vector(column) in the NQS
def np2nqs (npa,names=[]):
  nrow,ncol = numpy.shape(npa)
  if nrow < 2:
    print("np2nqs ERRA: must have at least 2 rows!")
    return None
  nqo = NQS(nrow)
  for i in range(nrow):
    vrow = py2vec( npa[i] )
    nqo.v[i].copy(vrow)
  for i in range(len(names)): # assign names
    nqo.s[i].s = names[i]
  return nqo

# converts nqs to numpy array
#  NB: each Vector(column) of nq is turned into a row in the numpy array returned
def nqs2np (nq,Full=False):
  if Full: nq.tog("DB")
  ncol,nrow = int(nq.m[0]),nq.size()
  npa = numpy.zeros( (nrow,ncol) )
  for i in range( ncol ):
    nptmp = vec2nparr( nq.getcol(nq.s[i].s) )
    npa[0:nrow][i] = nptmp[0:nrow]
  return npa

# convert nqs to a python dictionary
def nqs2pyd (nq,Full=False):
  d = {};
  if Full: nq.tog("DB");
  for i in range(int(nq.m[0])): d[nq.s[i].s] = vec2np(nq.getcol(nq.s[i].s))
  return d

# adds cols of nq2 to nq1
def NQAddToCols (nq1, nq2):
  cols = int(nq1.m[0])
  if cols != int(nq2.m[0]):
    print("unequal cols!")
    return False
  for i in range(cols):
    nq1.v[i].add(nq2.v[i])
  return True

# divide cols of nq by N
def NQDivCols (nq,N):
  cols = int(nq.m[0])
  for i in range(cols): nq.v[i].div(N)

