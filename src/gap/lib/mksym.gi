#mksym.gi
########################################################################
########################################################################
#
#
#
########################################################################
########################################################################

########################################################################
# 
#  Representations
# 
########################################################################
InstallGlobalFunction(Representations,
function(G,dim)
local k,i,j,tg,inds,realirrs,irrs,nofirrs,
parts,partition,thereps,res,combination,tmpres;
k:=0;res:=[];realirrs:=[];
tg:=CharacterTable(G);
irrs:=Irr(tg);
nofirrs:=Length(irrs);
inds:=Indicator(tg,irrs,2);



for i in [1..nofirrs]
do
if inds[i] = 1
then
realirrs[i]:=irrs[i];
elif inds[i]=0
then 
realirrs[i]:=irrs[i]+ComplexConjugate(irrs[i]);
elif inds[i]=-1
then
realirrs[i]:=2*irrs[i];
else
Print("ERROR!\n");
fi;
od;


parts:=Partitions(dim);
# parts:=[[dim]];
# only irreducible

for partition in parts
do
thereps:=List( partition, j -> 
Filtered(realirrs, chi -> (Degree(chi) = j) and (
Degree(chi) > DimFixedByGroup(chi)
 ) ));

	for combination in Cartesian(thereps)
	do
	tmpres:=Sum(combination);
	if Filtered( res, chi -> chi = tmpres) = []
	then
	k:=k+1;
	res[k]:=tmpres;
	fi;

	od;

od;
return res;
end );


#########################################################################
#
# HomomorphismsToSymmetric
#
#########################################################################

InstallGlobalFunction(HomomorphismsToSymmetric,
function(group,n)
# return the list of homomorphisms to the symmetric group
local FPiso, FPgroup, FPrelat, FPgens, gens, ns, N, p,Gbar,
is,iss,res,SG;
SG:=SymmetricGroup(n);

#FPiso:=IsomorphismFpGroup(group);
#FPgroup:=ImagesSource(FPiso);
#FPgens:=GeneratorsOfGroup(FPgroup);
#FPrelat:=RelatorsOfFpGroup(FPgroup);
# gens:=List(FPgens, t -> PreImagesRepresentative(FPiso,t)) ;                    

# begin temporary change

gens:=GeneratorsOfGroup(group);



ns:=Filtered(NormalSubgroups(group), N ->
        (Size(N) < Size(group) ) ); 
# nontrivial normal subgroups
is:=[];

for N in ns
do
p:=NaturalHomomorphismByNormalSubgroup(group,N);
Gbar:=ImagesSource(p);
iss:=IsomorphicSubgroups(SG,Gbar);
is[Position(ns,N)]:=List(iss, hom ->
 GroupHomomorphismByImages(group, SG, gens,
       List(gens, g ->
         Image(hom, Image(p, g) )
)));
 
od;                                                                            

# list of isomorphisms

res:=Flat(Filtered(is, t -> Length(t)>0 )); 
return res;

end  );




########################################################################
# 
#  MakeActions
# 
########################################################################
InstallGlobalFunction(MakeActions,
function(maxorbs,group,n,dim)
local i,j,k,Treps,nTreps,Vreps,nVreps, res, ii ,homs, nhoms,kersigma ,kerT,kerV,qq;

Treps:=Representations(group,2);
Vreps:=Representations(group,dim);
homs:=HomomorphismsToSymmetric(group,n);
nTreps:=Length(Treps);
nVreps:=Length(Vreps);                                                         
nhoms:=Length(homs);

ii:=0;
res:=[];

for i in [1..nTreps]
do
    # Print("Tchar number: ", i, "/", nTreps,"\n");
    kerT:=Kernel(Treps[i]);
    for j in [1..nVreps]
    do
        # Print("  - Vchar number: ", j, "/", nVreps,"\n");
	kerV:=Kernel(Vreps[j]);
	if (Size(Intersection(kerT,kerV))=1) 
	then

	for k in [1..nhoms]
	do
            # Print("    * nhom number: ", k, "/", nhoms,"\n");
	    kersigma:=Kernel(homs[k]);
		    qq:=LagSymmetryGroupCHARS(n,Treps[i],Vreps[j],homs[k]);
		    ## qq:=rec(n:=n,Tchar:=Treps[i],Vchar:=Vreps[j],sigma:=homs[k]);
            # Print(" [qq] ");

		if 	(Size(Intersection(kerT,kersigma)) = 1 ) and 
              		(Size(Intersection(kerV,kersigma)) <= 2 ) and  
			IsCoerciveCHARS(qq) # and 
			# (DimChiKerTCHARS(qq)>1) # and 
			and
			(not HasAlwaysCollisionsCHARS(qq)) 
		then
                    # Print(" [if] ");
		    ii:=ii+1;
		    res[ii]:=qq;
		    if ((maxorbs>0) and (ii=maxorbs))  then
		       return res;
		    fi;
		fi;
	    od;
	fi;
     od; 
od;
return res;
end);


########################################################################
# 
#  MakeBlockMatrix
# 
########################################################################


InstallGlobalFunction(MakeBlockMatrix,
function(matrix1,matrix2)
local i,j,result,l1,l2;
l1:=Length(matrix1);
l2:=Length(matrix2);
result:=List([1..(l1+l2)], a -> List([1..(l1+l2)], b -> 0 ));

for i in [1..l1] do
for j in [1..l1] do
		result[i][j] := matrix1[i][j];
od;
od;

for i in [1..l2] do
for j in [1..l2] do
		result[l1+i][l1+j] := matrix2[i][j];
od;
od;

return result;
end  ) ;


########################################################################
# 
#  MakeOneMatrix
# 
########################################################################


InstallGlobalFunction(MakeOneMatrix,
function(matrixlist)
local l,i,shortlist;
l:=Length(matrixlist);
if l=1 then
return matrixlist[1];
elif l=2 then
return MakeBlockMatrix(matrixlist[1],matrixlist[2]);
elif l>2 then
shortlist:=List([1..(l-1)], k -> matrixlist[k]);
return MakeBlockMatrix(MakeOneMatrix(shortlist),matrixlist[l]);
else
Error("MakeOneMatrix failed\n");
fi;
end ) ; 

########################################################################
# 
#  MatrixRepresentations
# 
########################################################################

InstallGlobalFunction(MatrixRepresentation,
function(Vchar)
local tg,irrs,irr_charvalues,irr_char,irr_hom,mult,i,
hom_constituents, hom_constituents_number,
gens,matrixgroup,listofmatrixes, g,irr_matrixes, group,
ind, already_there
;
group:=UnderlyingGroup(Vchar);
hom_constituents:=[];
hom_constituents_number:=0;
already_there:=false;

tg:=CharacterTable(group);
irrs:=IrreducibleRepresentationsByBaumClausen(group);
# irr_const:=ConstituentsOfCharacter(Vchar);

for irr_hom in irrs
do
	irr_charvalues:=List(ConjugacyClasses(group), 
		c -> TraceMat(Image(irr_hom,Representative(c))  ) );  
	irr_char:=Character(tg,irr_charvalues);
	ind:=Indicator(tg,[irr_char],2)[1];
	mult:=ScalarProduct(irr_char,Vchar);
	if (mult>0  and ind = 1) then
		for i in [1..mult] do
			hom_constituents_number:=hom_constituents_number + 1;
			if (Degree(irr_char) = 1) then
			hom_constituents[hom_constituents_number]:=irr_hom;
			elif (Degree(irr_char) = 2 ) then
			hom_constituents[hom_constituents_number]:=Realizza(group,irr_hom);
			else
			#hoping that it's with real coefficients!!
			hom_constituents[hom_constituents_number]:=Realizza(irr_hom);
			
			fi;
		od;
	elif (mult = 1 and ind = 0 ) then
	if (not already_there) then
			hom_constituents_number:=hom_constituents_number + 1;
			hom_constituents[hom_constituents_number]:=Raddoppia(group,irr_hom);
	# it needs to be of dimension 1... hence Raddoppia should work...
			already_there:=true;
	fi;
	elif (mult > 0 and ind <> 1) then
	Error ("I don't know what to do!\n mult > 1 and ind <> 1");
	# ind=1: done; ind=0: only mult=1; what about ind=-1?  not possible 
	# if there is one, it's of degree 1, hence not with ind=-1.
	# so everything else is a mistake.
	fi;

od;

gens:=GeneratorsOfGroup(group);

i:=0;
listofmatrixes:=[];
for g in gens do
i:=i+1;
irr_matrixes:=List( hom_constituents, phi -> Image (phi,g));
listofmatrixes[i]:=MakeOneMatrix(irr_matrixes);
od;

matrixgroup:=Group(listofmatrixes);

return GroupHomomorphismByImages(group, matrixgroup, gens, 
listofmatrixes);

end );


########################################################################
# 
#  Raddoppia
# 
########################################################################

InstallGlobalFunction(Raddoppia,
function(group,hom)
# WARNING: only in dim=1
local tg,gens,matrixes,matrixgroup,Xmatrix;
tg:=CharacterTable(group);
gens:=GeneratorsOfGroup(group);
Xmatrix := [ [ 1+E(4), -1 + E(4)], [1-E(4), -1-E(4) ] ];
matrixes:=List( gens, g -> 
( Xmatrix * 
MakeBlockMatrix(Image (hom,g), Image(hom,g^(-1)) )
* Xmatrix^(-1))
);

# don't need to check that it's working... 
# since it's called on 1-dimensional representations!!!
# if  Filtered( matrixes , g -> (not IsRealMatrix(g)) ) <> [] 
# then
# 	Error("Problem in Raddoppia!");
# 	return fail;
# fi;

matrixgroup:=Group(matrixes);

return GroupHomomorphismByImages(group, matrixgroup, gens, 
matrixes);

end )  ;



########################################################################
# 
#  MyDelta
# 
########################################################################

InstallGlobalFunction(MyDelta,
function(iii,jjj)
if iii=jjj then 
	return 1;
else 
	return 0;
fi; 
end )  ;


########################################################################
# 
#  IsUnitaryMatrix
# 
########################################################################
InstallGlobalFunction(IsUnitaryMatrix,
function(matrix)
local res, l, row, entry,i,j,k,tmpsum;
l:=Length(matrix);
# Print("*******\n\n");
for i in [1..l]
do
	for j in [1..l]
	do
	# Print ("i,j = ", i, ",", j , " ==> ");
		tmpsum:=0;
		for k in [1..l] do
		tmpsum := tmpsum +  matrix[i][k] * ComplexConjugate(matrix[j][k]);
		od; 
		# Print("v_i . v^*_j = ", tmpsum, "\n");
 		if tmpsum <> MyDelta(i,j)  then
			#Print("matrix not unitary with i,j =", i,j, "\n");
			#Print("tmpsum = ", tmpsum, "\n");
			#Print("MyDelta = ", MyDelta(i,j), "\n");
			return false;
		fi;
	od;
od;
return true;
end  );


########################################################################
# 
#  Realizza
# 
########################################################################
InstallGlobalFunction( Realizza,
function(group,hom)
local tg,gens,matrixes,matrixgroup,Xmatrix,listofmatrixes, tmpcheck;
### WARNING: ONLY WITH DIM=2 or if it's already real...!!!!
### what about stop checking that it's unitary? at the beginning only...

tg:=CharacterTable(group);
gens:=GeneratorsOfGroup(group);
Xmatrix := [ [ 1+E(4), -1 + E(4)], [1-E(4), -1-E(4) ] ];

listofmatrixes:=List( gens, g -> Image(hom,g) );

if Filtered( listofmatrixes, g -> not IsUnitaryMatrix(g) ) <>  []
then
	Print("filtered_listofmatrixes = ",
 Filtered( listofmatrixes, g -> not IsUnitaryMatrix(g) ) , "\n\n");

	Error("representation is not unitary!!!!");
	return fail;
fi;

matrixes:= Filtered( listofmatrixes , g -> not IsRealMatrix(g) );
if matrixes=[] then
	return hom;
else
	matrixes:=List( listofmatrixes, g -> 
		( Xmatrix * g * Xmatrix^(-1)));
	if  Filtered( matrixes , g -> not IsRealMatrix(g) ) <> [] 
	then
		Error("Problem in Realizza!");
		return fail;
	fi;

	matrixgroup:=Group(matrixes);

	return GroupHomomorphismByImages(group, matrixgroup, gens, 
	matrixes);
fi;
end ) ;



########################################################################
# 
#  IsReal
# 
########################################################################

InstallGlobalFunction(IsReal,
function(x)
return (ComplexConjugate(x) = x);
end ) ;

########################################################################
# 
#  IsRealMatrix
# 
########################################################################
InstallGlobalFunction(IsRealMatrix,
function(matrix)
local res, l, row, entry;
for row in matrix
do
	for entry in row
	do
		if not IsReal(entry) then
			return false;
		fi;
	od;
od;
return true;
end ) ;


########################################################################
# 
#  TypeOfAction
# 
########################################################################

InstallGlobalFunction( TypeOfAction,
function(qq)
local char,tmp;
char:=qq.Tchar;
tmp:=NumberOfReflections(char);
if (tmp = 0)
then
return 0;
elif (tmp = 1)
then
return 1;
elif (tmp>1)
then
return 2;
else
Error("NumberOfReflections not correct!\n");
return -1;
fi;
end ) ;


########################################################################
#
# ValueOfChar
#
########################################################################


InstallGlobalFunction( ValueOfChar,
function(char,elm)
local tg,ccs,pos,group;
tg:=UnderlyingCharacterTable(char);
group:=UnderlyingGroup(tg);
ccs:=ConjugacyClasses(tg);
pos:=Position(ccs,ConjugacyClass(group,elm));
return ValuesOfClassFunction(char)[pos];
end) ;

########################################################################
# 
#  MakeLSGfromCHARS
# 
########################################################################
InstallGlobalFunction(MakeLSGfromCHARS,
function(struct)
 local LSG,group,Tchar,Vchar,Vhom,sigma,n,tg,qq,action_type,
 dim,kernelTchar,p,Gbar,cyclic_order,cyclic_elms, tmpelms, cyclic_generator,
 kernel_gens,kern,el,perm,matrix,rotV,rotS,refV,refS;
 group:=struct!.group;
 n:=struct!.n;
 Tchar:=struct!.Tchar;
 Vchar:=struct!.Vchar;
 sigma:=struct!.sigma;
 Vhom:=MatrixRepresentation(Vchar);
 tg:=CharacterTable(group);
 qq:=rec(n:=n,Tchar:=Tchar,Vchar:=Vchar,sigma:=sigma);
 action_type:=TypeOfAction(qq);
 dim:=struct!.dim;
 kernelTchar:=Kernel(Tchar);
 p:=NaturalHomomorphismByNormalSubgroup(group,kernelTchar);
 Gbar:=ImagesSource(p);
 if IsCyclic(Gbar)
  then
  # in case Gbar cyclic
   cyclic_order:=Size(Gbar);
   cyclic_elms:=List(Elements(Gbar), gb ->  PreImagesRepresentative(p,gb));
   tmpelms:=Filtered(cyclic_elms, gg -> Order(Image(p,gg))=cyclic_order);
  else
   cyclic_order:=Size(Gbar)/2;
   cyclic_elms:=Filtered(List(Elements(Gbar), gb ->  PreImagesRepresentative(p,gb)), gg -> ValueOfChar(DeterminantOfCharacter(Tchar),gg)=1);
# return the elements of Gbar (actually their preimages in G)
# that act orientation-preserving on T
   tmpelms:=Filtered(cyclic_elms, gg -> Order(Image(p,gg))=cyclic_order);
fi;

## 
if (tmpelms <> [])  then
cyclic_generator:=tmpelms[1];
else
Print("tmpelms = ", tmpelms, "\n");
Error("Problem: cyclic_generator not found.\n");
return;
fi;
## TODO!!!
 kernelTchar:=Kernel(Tchar);
 kernel_gens:=[];
 for el in GeneratorsOfGroup(kernelTchar)
   do
     perm:=Image(sigma,el);
     matrix:=Image(Vhom,el);
     Append(kernel_gens, [Tuple([ matrix,perm ] )] );
   od;
 kern:=GroupWithGenerators(kernel_gens);  
## enter kernel...

rotV:= Image(Vhom,cyclic_generator);
rotS:= Image(sigma,cyclic_generator);

refV:=[];
refS:=();

if (action_type > 0 ) then
  tmpelms:=Filtered(List(Elements(Gbar), gb ->  PreImagesRepresentative(p,gb)), gg -> ValueOfChar(DeterminantOfCharacter(Tchar),gg)=-1);


 if (Length(tmpelms) <= 1 ) then
   Error("Problem: 2 non_cyclic elements not found.\n");
   return;
 fi;

refV:=Image(Vhom,tmpelms[1]);
refS:=Image(sigma,tmpelms[1]);

fi;

## now build LagSymmetryGroup...
# action_type, n, kern, rotV, rotS, 

LSG:=LagSymmetryGroup(action_type,n , kern, rotV, rotS, refV, refS);

return LSG;
end
);




########################################################################
# 
#  MkOrbs
# 
########################################################################
InstallGlobalFunction( MkOrbs,
function(maxorbs,dirname, order, n,dim)
local ind,group,aa,qq,filename,stru,tmporb,LSG;
ind:=0;

for group in AllGroups(Size,order) do
tmporb:=0;
        Print("IDGroup = ", IdGroup(group), "\n");
        aa:=MakeActions(maxorbs,group,n,dim);
        for qq in aa do
        tmporb:=tmporb+1;
        ind:=ind+1;
        Print("ind=", ind, "; tmporb=", tmporb, ";\n");
        filename:=Concatenation(dirname, "/orbit_");
        Append(filename,String(n));
        Append(filename,"b_");
        Append(filename,String(order));
        Append(filename,"g_");
        Append(filename,String(ind));
        Print("generating LSG...\n");
	LSG:=MakeLSGfromCHARS(qq);
        Print("creating file...", filename , ";\n");
        MakeMinorbSymFile(filename,LSG);
        Print("...done!\n");
        od;
od;
UnloadSmallGroupsData();
end );
