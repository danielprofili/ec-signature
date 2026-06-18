with(LinearAlgebra):
with(combinat):
with(Iterator):
with(ListTools):

svp := proc(f,v)   
# Sign Variation of Polynomial
# In : f, a polynomial in v with real coeffcients
#      v, a variable
# Out: the sign variataion count of the coefficients of f
  local d,L,i;
  d  := degree(f,v);
  L  := [coeff(f,v,d-i)$i=0..d];
  return svl(L);
end:


svl := proc(L)    
# Sign Variation of List
# In : L, a list of real numbers
# Out: the sign variation count on L
  local Lb, Cb, i, c;
  Lb := remove(v->evalf(v=0), L);
  Cb := [seq(Lb[i]*Lb[i+1],i=1..nops(Lb)-1)];
  Cb := remove(v->evalf(v>0), Cb);
  c  := nops(Cb);
  return c;
end:  


cpl := proc(Ls)    
# Cartesian Product of Lists
# In : Ls, a list of lists
# Out: The cartesian product of Ls
  local P,C;
  P := cartprod(Ls);
  C := [];
  while not P[finished] do
    C := [op(C), P[nextvalue]()];
  end do;
  
  # C := [seq(Reverse(convert(p[],list)),p=P)];
  return C;
end:


kpm := proc(Ms)    
# Kroneker Product of Matrices
# In : Ms, a list of matrices
# Out: The Kronecker products of Ms
  return foldl(KroneckerProduct,op(Ms));
end:


kpl := proc(Ls)    
# Kronecker Product of Lists
# In : Ls, a list of lists
# Out: The Kronecker products of Ls
  local Ms,L;
  Ms := map(L-><op(L)>,Ls);
  return map(expand,convert(kpm(Ms),list));
end:


sig := proc(M)      
# Signature of Matrix 
# In : M, a real symmetric matrix
# Out: the signature (difference between number of positive and negative eigenvalues) of M
  local cp,n,s;
  n  := RowDimension(M);
  cp := CharacteristicPolynomial(M,lambda);
  s  := 2*svp(cp,lambda) - n + ldegree(cp,lambda);
  return s;
end:


epm := proc(f,v,M)   
# Evaluate Polynomial on Matrix
# In : f, a polynomial in v
#      v, a variable
#      M, a matrix
# Out: f(M)
  local F,fe,k;
  fe := expand(f);
  F  := simplify(coeff(fe,v,0)*IdentityMatrix(RowDimension(M)) + add(coeff(fe,v,k)*M^k,k=1..degree(fe)));
  return F;
end:

esvp := proc(f,v)   
# Extended Sign Variation of Polynomial
# In : f, a polynomial in v with real coefficients
#      v, a variable
# Out: the sign variataion count on the coefficients of f
  local d,L,i,c;
  d  := degree(f,v);
  L  := [coeff(f,v,d-i)$i=0..d];
  c  := esvl(L);
  return c;
end:

esvl := proc(L)  
# Extended Sign Variation of List
# In : L, a list of real numbers
# Out: extended sign variation count on L
  local inds,es,svcs,Le,c,i,e;
  inds := [SearchAll(0, L)];
  es := cpl([[-1,1]$nops(inds)]);
  svcs := 0;
  for e in es do
    Le := L;
    for i from 1 to nops(inds) do
      Le[inds[i]] := e[i];
    od:
    svcs := svcs + svl(Le);
  od:
  
  c := 1/nops(es) * svcs;
  return c; 
end:

V := proc(m)      
# V matrix
# In :  m, positive integer
# Out:  V, m by 2^m V matrix
  local Vm;
  Vm  := cpl([[-1,1]$m]);
  Vm  := map(v->svl([op(v),1]),Vm);
  Vm  := Matrix(m,2^m,(t,s)->`if`(Vm[s]=m-t,1,0));
  return Vm;
end:

H := proc(m)     
# H matrix
# In : m, positive integer
# Out: H, a 2^m by 2^m H matrix
  return kpm([Matrix(2,2,[[1,1],[-1,1]])$m]);
end:

Csig := proc(m)   
# Combinatorial part for signature method
# In : m, positive integer
# Out: Csig
  return V(m) . H(m)^(-1);
end:

EC := proc(F,G)    
# Eigenvalue Configuration. signature method
# In : F, G, numeric real symmetric matrices
# Out: Eigenvalue configuration of F and G   
  local m,n,f,fp,es,fes,i,feGs,hes,sigma,Vm,Hm,c,e,fe,feG,sves,he,tdes,r;
  m     := RowDimension(F);
  n     := RowDimension(G);  
  f     := CharacteristicPolynomial(F,x);
  fp    := [seq(diff(f,[x$i]),i=0..m-1)];
  es    := cpl([[0,1]$m]);
  fes   := [seq(mul(fp[i+1]^e[i+1],i=0..m-1),e=es)];
  feGs  := [seq(epm(fe,x,G),fe=fes)];
  hes   := [seq(CharacteristicPolynomial(feG,x),feG=feGs)];
  sves  := [seq(svp(he,x),he=hes)];
  tdes  := [seq(ldegree(he,x),he=hes)];
  sigma := [seq(2*sves[i]-n+tdes[i],i=1..nops(sves))];
  sigma := Vector(sigma);
  Hm    := H(m);
  Vm    := V(m);
  r     := Hm^(-1).sigma;
  c     := Vm.r;
  return c;
end:

Asig := proc(F,G)   
# Algebraic part for the signature method
# In : F, G symbolic real symmetric matrices
# Out: Asig(F,G)
  local m,f,fp,i,es,fes,e,feGs,fe,hes,feG;
  m := RowDimension(F);
  #F := Matrix(m,m,(i,j)->a[i,j], shape=symmetric);
  #G := Matrix(n,n,(i,j)->b[i,j], shape=symmetric);
  f := CharacteristicPolynomial(F, x);
  
  fp := [seq(diff(f,[x$i]),i=0..m-1)];
  es    := cpl([[0,1]$m]);
  fes   := [seq(mul(fp[i+1]^e[i+1],i=0..m-1),e=es)];
  feGs  := [seq(epm(fe,x,G),fe=fes)];
  hes   := [seq(CharacteristicPolynomial(feG,x),feG=feGs)];
  return Vector(hes);
end:

ECSig := proc(F,G)   
# Conditions for EC, signature method
# In : F, G, symbolic real symmetric matrices
# Out: Csig, Asig such that c = EC(F,G) iff c = Csig . sig(Asig)
  return Csig(RowDimension(F)), Asig(F,G);
end:
