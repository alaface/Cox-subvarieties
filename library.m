////////////////////////////////////////////////
// NewGens
// INPUT = Ideal of subvariety, Basis of \tilde I_1, grading map
// OUTPUT = new generators 
////////////////////////////////////////////////

NewGens := function(I,lis,h)
 R1 := Parent(lis[1]);
 R2 := Parent(I.1);
 r1 := Rank(R1);
 lis := [Evaluate(f,[R2.i : i in [1..r1]]) : f in lis];
 den := lis[1];
 n := 0;
 repeat
  n := n + 1;
  ide := I+Ideal(den)^n;
  J := &meet[Saturation(ide,I+Ideal(b)) : b in lis[2..#lis]];
  ll := [f : f in MinimalBasis(J) | f notin ide];
 until #ll gt 0;
 E1 := Domain(h);
 Cl := Codomain(h);
 l := [h(b) : b in Basis(E1)];
 E2 := RSpace(BaseRing(E1),Rank(E1) + #ll);
 w := [h(E1!Exponents(Monomials(f)[1]) - n*E1!Exponents(den)) : f in ll];
 h := hom<E2->Cl|l cat w>;
 return ll,lis,n,h;
end function;


////////////////////////////////////////////////
// NewRing
// INPUT = Ideal of subvariety, ll, m, n
// OUTPUT = the ideal of the new ring
////////////////////////////////////////////////

NewRing := function(I,ll,lis,n,h)
 den := lis[1];
 R1 := Parent(ll[1]);
 Q := BaseRing(R1);
 r1 := Rank(R1);
 r2 := r1 + #ll;
 R2 := PolynomialRing(Q,r2);
 v := [R2.i : i in [1..r1]];
 den := Evaluate(den,v);
 ll := [Evaluate(f,v) : f in ll];
 lis := [Evaluate(f,v) : f in lis];
 I := Ideal([Evaluate(g,v) : g in MinimalBasis(I)]);
 I := Saturation(I + Ideal([R2.(r1+i)*den^n-ll[i] : i in [1..#ll]]),den);
 ln := [];
 for f in Basis(I) do
  l := [m : m in Monomials(f) | Degree(m) eq 1];
  if #l gt 0 then Append(~ln,l[1]); 
  end if;
 end for;
 ln := Setseq(Set(ln));
 ind := [i : i in [1..r2] | R2.i notin ln];
 I1 := EliminationIdeal(I,{R2.i : i in ind});
 I2 := EliminationIdeal(I+Ideal(lis),{R2.i : i in ind});
 R3<[T]> := PolynomialRing(Q,r2 - #ln);
 u := [R3!1 : i in [1..r2]];
 E1 := Domain(h);
 Cl := Codomain(h);
 E2 := RSpace(BaseRing(E1),#ind);
 h := hom<E2->Cl|[h(E1.i) : i in ind]>;
 for i in [1..r2] do
  if i in ind then 
   j := Position(ind,i);
   u[i] := T[j];
  end if;
 end for;
 I1 := Ideal([Evaluate(f,u) : f in Basis(I1)]);
 I2 := Ideal([Evaluate(f,u) : f in Basis(I2)]);
 return I1,MinimalBasis(I2),h;
end function;


////////////////////////////////////////////////
// IsCox
// INPUT = Ideal of subvariety, Basis of \tilde I_1
// OUTPUT = boolean
////////////////////////////////////////////////

IsCox := function(I,lis)
 R1 := Parent(lis[1]);
 R2 := Parent(I.1);
 r1 := Rank(R1);
 r2 := Rank(R2);
 lis := [Evaluate(f,[R2.i : i in [1..r1]]) : f in lis];
 codim := Dimension(I) - Dimension(I+Ideal(lis));
 return codim gt 1;
end function;


////////////////////////////////////////////////
// ComputeCox
// INPUT = Ideal of subvariety, Basis of \tilde I_1, grading map
// OUTPUT = boolean
////////////////////////////////////////////////

ComputeCox := function(I,lis,h)
 repeat
  ll,lis,n,h := NewGens(I,lis,h);
  I,lis,h := NewRing(I,ll,lis,n,h);
 until IsCox(I,lis);
 gradings := Transpose(Matrix([h(b) : b in Basis(Domain(h))]));
 return MinimalBasis(I),gradings;
end function;


////////////////////////////////////////////////
// IdeCodim2
// INPUT = Toric variety
// OUTPUT = codimension 2 irrelevant locus
////////////////////////////////////////////////

IdeCodim2 := function(Z)
 I := IrrelevantIdeal(Z);
 r := Rank(Parent(I.1));
 return Basis(&meet[J : J in RadicalDecomposition(I) | r-Dimension(J) eq 2]);
end function;

