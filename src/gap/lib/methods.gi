InstallGlobalFunction(cos2piover,
function(n)
return (E(n) + E(n)^(-1))/2;
end);

InstallGlobalFunction(sin2piover,
function(n) return (-E(4))*(E(n) - E(n)^(-1))/2;
end);

InstallGlobalFunction(RotationMatrixDim2,
function(order)
return
[ [ cos2piover(order), -sin2piover(order)],
  [ sin2piover(order), cos2piover(order) ]
   ];
end);


InstallGlobalFunction(RotationMatrixDim3,
function(order)
return
[ [ cos2piover(order), -sin2piover(order),0],
  [ sin2piover(order), cos2piover(order),0 ],
  [ 0 , 0, 1] ];
end);

###########################
###########################
###########################
###########################

InstallMethod(SRC_QQ,"for LSGs" , [IsLagSymmetryGroup],
function(LSG)
local gens,group,tg,Vchar,Tchar,sigma,listofVtraces,x,
charvalues,c,Tgens,cyclic_order,glob_gens,i;


cyclic_order:= Size ( Group( Concatenation(GeneratorsOfGroup(LSG!.kerneltau),[LSG!.rotation]) )) / Size(LSG!.kerneltau);

Tgens:= List( GeneratorsOfGroup(LSG!.kerneltau) , x -> DiagonalMat([1,1]) );
Append(Tgens,[ RotationMatrixDim2(cyclic_order) ] );

if LSG!.action_type=0 then
gens:=Concatenation(
    GeneratorsOfGroup(LSG!.kerneltau),
	[LSG!.rotation]
	);
else
gens:=Concatenation(
    GeneratorsOfGroup(LSG!.kerneltau),
	[LSG!.rotation,LSG!.reflection]
	);
    Append(Tgens, [  [ [0,1],[1,0] ] ]);
fi;


glob_gens:=List(
[1..Length(Tgens)], i->
Tuple( [ Tgens[i], gens[i][1], gens[i][2] ] )
);

group:=GroupWithGenerators( glob_gens );
tg:=CharacterTable(group);

sigma := GroupHomomorphismByImages( group, SymmetricGroup( LSG!.NOB ),
	glob_gens,
		List(gens, x -> x[2])
	);

Vchar:=Character(tg,
	List(  ConjugacyClasses( group ), function ( c )
            return TraceMat(  Representative(c)[2] ); end )
);



Tchar:= Character(tg,
	List(  ConjugacyClasses( group ), function ( c )
            return TraceMat( Representative(c)[1]  ); end )
);

## return LagSymmetryGroupCHARS( LSG!.NOB, Tchar, Vchar, sigma );

##!
return rec(
	Vchar:=Vchar,
	Tchar:=Tchar,
	sigma:=sigma,
	n:=LSG!.NOB,
	group:=group
	  );
end);




###
InstallGlobalFunction(
ConfigurationSpaceChar,
 function ( group, n, dim, Vchar, sigma )
    local  tg, characters, phi, SG, k, single, tensored_char, sigmachar;
    tg := CharacterTable( group );
    SG := SymmetricGroup( n );
    sigmachar:= RestrictedClassFunctions( [ NaturalCharacter( SG ) - TrivialCharacter( SG ) ], sigma )[1] ;
    single := [];
    for k in [ 1 .. Length( Vchar )] do
    single[k] := Vchar[k] * sigmachar[k];
    od;
    tensored_char:= VirtualCharacter(tg,single);
    return tensored_char;
end);

##
InstallGlobalFunction(DimFixed,
function ( tg, H, chi )
    local  th;
    if IsTrivial( H )  then
        return Degree( chi );
    else
        th := CharacterTable( H );
        FusionConjugacyClasses( th, tg );
        return ScalarProduct( th, TrivialCharacter( th ), RestrictedClassFunction( tg, chi, H ) );
    fi;
    return;
end);


###
InstallGlobalFunction(
DimFixedByGroup, function ( char )
    return DimFixed( UnderlyingCharacterTable( char ), UnderlyingGroup( char ), char );
end);


##
InstallGlobalFunction(DimChiKerT,
 function ( LSG )
    local  qq,group, n, dim, Vchar, sigma, kerT, char;
    qq:=Objectify(LagSymmetryGroupCHARSType,SRC_QQ(LSG));
    return DimChiKerTCHARS(qq);


    Vchar := qq.Vchar;
    sigma := qq.sigma;
    dim := Degree( Vchar );
    n := qq.n;
    kerT := Kernel( qq.Tchar );
    group := UnderlyingGroup( Vchar );
    char := ConfigurationSpaceChar( group, n, dim, Vchar, sigma );
    return DimFixed( UnderlyingCharacterTable( char ), kerT, char );

 end);



##

########################################################################
InstallGlobalFunction(DimChiKerTCHARS,
 function ( qq )
    local  group, n, dim, Vchar, sigma, kerT, char;
    Vchar := qq!.Vchar;
    sigma := qq!.sigma;
    dim := Degree( Vchar );
    n := qq!.n;
    kerT := Kernel( qq!.Tchar );
    group := UnderlyingGroup( Vchar );
    char := ConfigurationSpaceChar( group, n, dim, Vchar, sigma );
    return DimFixed( UnderlyingCharacterTable( char ), kerT, char );
 end);


#




###
InstallGlobalFunction( IsotropySubgroups,
function ( G, char )
    local  tm, AllTheSubgroups, l, TOMMatrix, tg, dimsfixed, Isotropies;
	if Size(G) = 1 then
		return [];
	fi;
    tg := CharacterTable( G );
    tm := TableOfMarks( G );
    TOMMatrix := MatTom( tm );
    l := Length( SubsTom( tm ) );
    AllTheSubgroups := List( [ 1 .. l ], function ( k )
            return RepresentativeTom( tm, k );
        end );
    tg := CharacterTable( G );
    dimsfixed := List( AllTheSubgroups, function ( H )
            return DimFixed( tg, H, char );
        end );
    Isotropies := Filtered( [ 1 .. l ], function ( j )
            return Length( Filtered( [ 1 .. l ], function ( i )
                        return (TOMMatrix[i][j] > 0 and dimsfixed[i] = dimsfixed[j] and dimsfixed[j] > 0);
                    end ) ) = 1;
        end );
    return List( Isotropies, function ( h )
            return AllTheSubgroups[h];
        end );
end);



InstallGlobalFunction( MaximalIsotropies ,
function ( group, isotropy_subgroups )
    local  i, K, result, H, is_maximal;
    result := [  ];
    i := 1;
    for H  in isotropy_subgroups  do
        is_maximal := true;
        for K  in isotropy_subgroups  do
            if not H = K and not Filtered( ConjugateSubgroups( group, K ), function ( k )
                           return IsSubgroup( k, H );
                       end ) = [  ]  then
                is_maximal := false;
                break;
            fi;
        od;
        if is_maximal  then
            result[i] := H;
            i := i + 1;
        fi;
    od;
    return result;
end);


###
InstallGlobalFunction( NumberOfReflections,
 function ( char )
    local  res, h, isotropies, maximal_isotropies, group;
    group := UnderlyingGroup( char );
    isotropies := IsotropySubgroups( group, char );
    maximal_isotropies := MaximalIsotropies( group, isotropies );
    res := -1;
    if Length( isotropies ) = 1  then
        res := 0;
    elif Length( isotropies ) >= 2  then
        res := Sum( List( maximal_isotropies, function ( h )
                  return Length( ConjugateSubgroups( group, h ) );
              end ) );
    else
        res:=0;
        ## Error( " Length(isotropies)= ", Length( isotropies ), "\n!!!" );
    fi;
    return res;
end);



########
InstallGlobalFunction( IsTypeRDirection,
function (LSG,dir)
local all_matrixes,small_matrixes,row,column,x,matrix,qq,group,glob_gens,
tg,tmpVhom,tmpVchar,c,matrix_group;
if ( (not dir in [1,2,3] ) or (not LSG!.dim = 3) )
then
	Print("Error: IsTypeRDirection wrong argument\n");
	return fail;
fi;


all_matrixes:=List(GeneratorsOfGroup(LSG!.kerneltau) , x -> x[1]);
Append(all_matrixes, [ LSG!.rotation[1] ]);
if LSG!.action_type>0
then
	Append(all_matrixes, [ LSG!.reflection[1] ]);
fi;

if not ForAll( all_matrixes,
# test function
    function (matrix)
    local row, column;
    row:=matrix[dir];
    column:=List(matrix, x -> x[dir]);
    return ( Filtered( [1,2,3], x -> x=dir or row[x]=0) = [1,2,3]  ) and
    ( Filtered( [1,2,3], x -> x=dir or row[x]=0) = [1,2,3]  );
    end )
then
	return false;
fi;
# define the representation
#
# cut the dir direction first
small_matrixes:=List( all_matrixes ,   function(mat)
local tmp, tmp_row, tmp_column, i, j ;
tmp:= [  ];
for i in Filtered( [1,2,3] , x -> x<> dir)
do

	tmp_row:=[];
	for j in Filtered( [1,2,3] , y ->  y <> dir)
	do
		Append(tmp_row, [ mat[i][j] ] );
	od;
	Append(tmp,[ tmp_row ] );
od;
return tmp;
end);

qq:=SRC_QQ(LSG);
group:=qq.group;
tg:=CharacterTable(group);
glob_gens:=GeneratorsOfGroup(group);
matrix_group:=GroupWithGenerators(small_matrixes);

tmpVhom:=GroupHomomorphismByImages( group, matrix_group,
        glob_gens, small_matrixes );

tmpVchar:=Character(tg,
        List(  ConjugacyClasses( group ), function ( c )
            return TraceMat(  Image(tmpVhom,Representative(c)) ); end )
);
return tmpVchar = qq.Tchar;

end);


########
InstallGlobalFunction( TypeRDirections,
function(LSG)
if  (LSG!.dim = 2)
then
    if IsTypeR(LSG) then
	    return [3];
    else
	    return [];
    fi;
elif ( LSG!.dim = 3 )
then
	return Filtered( [1,2,3], x -> IsTypeRDirection(LSG,x) );
else
	return fail;
fi;
end);

####

InstallGlobalFunction(TransitiveDecomposition,
function(LSG)
local sims;
sims:=List(GeneratorsOfGroup(LSG!.kerneltau) , h -> h[2] );
Append(sims, [ LSG!.rotation[2] ]);
if LSG!.action_type > 0
then
Append(sims, [LSG!.reflection[2] ]);
fi;
return Orbits( Group ( sims ) , [1..LSG!.NOB]);
end);






########

InstallMethod(IsValidLSG,"for LSGs",[IsLagSymmetryGroup],
function(LSG)
local gens,group,kertau,proj,Gbar,Rotbar,Refbar;
if LSG!.action_type=0 then
    gens:=Concatenation(
	GeneratorsOfGroup(LSG!.kerneltau),
	    [LSG!.rotation]
	    );
else
    gens:=Concatenation(
	GeneratorsOfGroup(LSG!.kerneltau),
	    [LSG!.rotation,LSG!.reflection]
	    );
fi;

group:=GroupWithGenerators( gens );
kertau:=Subgroup(group,GeneratorsOfGroup(LSG!.kerneltau));
## check that kertau is normal.
if ( not IsNormal(group,kertau)  )
then
	return false;
elif ( LSG!.action_type = 0 ) then
	return true;
fi;
## so far, if action_type=0 it's done...


## now we assume it is normal
if Size(kertau)>1 then
    proj:=NaturalHomomorphismByNormalSubgroup(group,kertau);
    Gbar:=ImagesSource(proj);
    Rotbar:=Image(proj,LSG!.rotation);
    Refbar:=Image(proj,LSG!.reflection);
else
    Gbar:=group;
    Rotbar:=LSG!.rotation;
    Refbar:=LSG!.reflection;
fi;


if ( 	( Refbar^2 =  Identity(Gbar) )
	and
	( Refbar^(-1) * Rotbar * Refbar = Rotbar^(-1) )
   )
then
	return true;
else
	return false;
fi;

end);


## true methods


InstallMethod(IsCoercive,"for LSGs",[IsLagSymmetryGroup],
function(LSG)
    local qq,group,n,dim,Vchar,sigma;
    qq:=SRC_QQ(LSG);
    Vchar:=qq.Vchar;
    sigma:=qq.sigma;
    dim:=Degree(Vchar);
    n:=LSG!.NOB;
    group:=UnderlyingGroup(Vchar);
    return DimFixedByGroup(
    ConfigurationSpaceChar(group,n,dim,Vchar,sigma)
    ) = 0;
end);

InstallMethod(IsCoerciveCHARS,"for LSG-CHARS",[IsLagSymmetryGroupCHARS],
function(qq)
    local group,n,dim,Vchar,sigma;
    Vchar:=qq!.Vchar;
    sigma:=qq!.sigma;
    dim:=Degree(Vchar);
    n:=qq!.n;
    group:=UnderlyingGroup(Vchar);
    return DimFixedByGroup(
    ConfigurationSpaceChar(group,n,dim,Vchar,sigma)
    ) = 0;
end);


InstallMethod(HasAlwaysCollisions, "for LSGs", [IsLagSymmetryGroup],
function ( LSG )
local qq;
qq:=Objectify(LagSymmetryGroupCHARSType, SRC_QQ(LSG) ) ;
return HasAlwaysCollisionsCHARS(qq);
end);



InstallMethod(HasAlwaysCollisionsCHARS, "for LSG-CHARSs", [IsLagSymmetryGroupCHARS],
function(qq)
local Tchar,Vchar,hom,maximal_isotropies,group,isotropies,
n,tg,gens,group_times_SG,tgts,g,myPhi,Vrep_times,chira,
is_only_collision,
H,h,SG,Sigma_rep_times, i,j,k,
gens_of_H,
Phi_of_gens_of_H,
Phi_of_H,Phi_of_H_tilde,
ll,dim1,dim2,tmplist,Ogroup,Vhom;
# Print("Beginning HasAlwaysCollisionsCHARS... \n");

Tchar:=qq!.Tchar;
Vchar:=qq!.Vchar;
hom:=qq!.sigma;
n:=qq!.n;
group:=UnderlyingGroup(Tchar);

tg:=CharacterTable(group);
gens:=GeneratorsOfGroup(group);

Vhom:=MatrixRepresentation(Vchar);
Ogroup:=Image(Vhom);
SG:=SymmetricGroup(n);

isotropies:=IsotropySubgroups(group,Tchar);
maximal_isotropies:=MaximalIsotropies(group,isotropies);
# Print("Maximal Isotropies found! \n");

group_times_SG:=DirectProduct(Ogroup,SG);
myPhi:=GroupHomomorphismByImages( group, group_times_SG, gens,
       List( gens, function(g)
                       return   Tuple( [ Image(Vhom,g) , Image( hom, g) ] )  ;
		  end
           )
	                          );



# Vrep_times:=Character(CharacterTable(Ogroup),
#         List(  ConjugacyClasses( Ogroup ), function ( c )
#	            return TraceMat(  Representative(c) ); end )
#		    );

# Vrep_times:=RestrictedClassFunction( Vchar, Vhom);

# Vrep_times:=RestrictedClassFunction( Vchar,  Projection(group_times_SG,1) );

# Sigma_rep_times:= NaturalCharacter(SG);
## RestrictedClassFunction(NaturalCharacter(SG),Projection(group_times_SG,2));

# chira:=Tensored([Vrep_times],[Sigma_rep_times])[1];


# scan the ordered pairs i<j
is_only_collision:=false;
for i in [1..(n-1)] do
# Print("i=",i,"; " );
if is_only_collision then break; fi;
	for j in [(i+1)..n] do
        # Print("j=",j,"\n");
if is_only_collision then break; fi;

for H in maximal_isotropies do
gens_of_H:=GeneratorsOfGroup(H);
Phi_of_gens_of_H:=List( gens_of_H, h-> Image(myPhi,h) );
# Phi_of_H:= Subgroup(group_times_SG,Phi_of_gens_of_H);
Phi_of_H:=Image(myPhi,H);

ll:=Phi_of_gens_of_H;

Append(ll, [Image(Embedding(group_times_SG,2),(i,j) ) ]);

Phi_of_H_tilde:=Subgroup(group_times_SG, ll);

# dim1:=DimFixed(tgts,Phi_of_H,chira);
# dim2:=DimFixed(tgts,Phi_of_H_tilde,chira);

dim1:= Sum( Elements(Phi_of_H) , x -> TraceMat(  x[1] ) * (n - NrMovedPoints( x[2] )  ) )   / Size(Phi_of_H) ;
dim2:= Sum( Elements(Phi_of_H_tilde) , x -> TraceMat(  x[1] ) * (n - NrMovedPoints( x[2] )  ) )   / Size(Phi_of_H_tilde) ;

if (dim1 = dim2 ) then
is_only_collision:=true;
break;
fi;
od; #end of for H

od; # end for j
od; # end for i

return is_only_collision;
end);




InstallMethod(HasAlwaysCollisionsCHARStemp, "for LSG-CHARSs", [IsLagSymmetryGroupCHARS],
function ( qq )
local  Tchar, Vchar, hom, maximal_isotropies, group, isotropies, n,
tg, gens, group_times_SG, tgts, g, myPhi, Vrep_times,
chira, is_only_collision, H, h, SG, Sigma_rep_times,
i, j, k, gens_of_H, Phi_of_gens_of_H, Phi_of_H, Phi_of_H_tilde, ll, dim1,
dim2;
## qq:=SRC_QQ(LSG);
    Tchar := qq!.Tchar;
    Vchar := qq!.Vchar;
    hom := qq!.sigma;
    n := qq!.n;
    group := UnderlyingGroup( Tchar );
    SG := SymmetricGroup( n );
    isotropies := IsotropySubgroups( group, Tchar );
    maximal_isotropies := MaximalIsotropies( group, isotropies );
    tg := CharacterTable( group );
    gens := GeneratorsOfGroup( group );
##     group_times_SG := DirectProduct( group, SG );
##     tgts := CharacterTable( group_times_SG );
    group_times_SG:=group;
    tgts:=tg;
##    myPhi := GroupHomomorphismByImages( group, group_times_SG, gens, List( gens, function ( g )
##              return Image( Embedding( group_times_SG, 1 ), g ) * Image( Embedding( group_times_SG, 2 ), Image( hom, g ) );
##          end ) );
      myPhi:=GroupHomomorphismByImages( group, group, gens,
       List( gens, function(g)
                       return   Tuple( [ g[1], Image( hom, g) ] )  ;
                       #return Tuple([ g[1] ,Image( hom , g ) ] );
		  end
           )
	                          );

##     Vrep_times := RestrictedClassFunction( Vchar, Projection( group_times_SG, 1 ) );
    Vrep_times := Vchar;
    Sigma_rep_times := RestrictedClassFunction( NaturalCharacter( SG ), Projection( group_times_SG, 2 ) );
    chira := Tensored( [ Vrep_times ], [ Sigma_rep_times ] )[1];
    is_only_collision := false;
    for i  in [ 1 .. n - 1 ]  do
        if is_only_collision  then
            break;
        fi;
        for j  in [ i + 1 .. n ]  do
            if is_only_collision  then
                break;
            fi;
            for H  in maximal_isotropies  do
                gens_of_H := GeneratorsOfGroup( H );
                Phi_of_gens_of_H := List( gens_of_H, function ( h )
                        return Image( myPhi, h );
                    end );
                Phi_of_H := Subgroup( group_times_SG, Phi_of_gens_of_H );
                ll := Phi_of_gens_of_H;
                Append( ll, [ Image( Embedding( group_times_SG, 2 ), (i,j) ) ] );
                Phi_of_H_tilde := Subgroup( group_times_SG, ll );
                dim1 := DimFixed( tgts, Phi_of_H, chira );
                dim2 := DimFixed( tgts, Phi_of_H_tilde, chira );
                if dim1 = dim2  then
                    is_only_collision := true;
                    break;
                fi;
            od;
        od;
    od;
    return is_only_collision;
end);


InstallMethod(IsTypeR,"for LSGs",[IsLagSymmetryGroup],
function(LSG)
    local qq,group,detVchar,detTchar;
    if ( LSG!.dim = 2 )
    then
	qq:=SRC_QQ(LSG);
	detVchar:=DeterminantOfCharacter(qq.Vchar);
	detTchar:=DeterminantOfCharacter(qq.Tchar);
    	return detVchar = detTchar;
    elif ( LSG!.dim = 3 )
    then
	return Filtered( [1,2,3], x -> IsTypeRDirection(LSG,x) ) <> [];
    else
		return fail;
    fi;
end);


InstallMethod(IsTransitiveLSG,"for LSGs",[IsLagSymmetryGroup],
function(LSG)
    return Length(TransitiveDecomposition(LSG)) = 1;
end);


InstallMethod(ActionType,"for LSGs",[IsLagSymmetryGroup],
function(LSG)
local char,tmp;
 char:=SRC_QQ(LSG).Tchar;
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
 return fail;
 fi;
end);



##
# IsRedundant:=NewProperty("IsRedundant",IsLagSymmetryGroup);

InstallMethod(IsRedundant,"for LSGs", [IsLagSymmetryGroup],
function(LSG)
  local qq;
  if LSG!.dim=2
  then
    if ( not IsTypeR(LSG) )
    then
      qq:=SRC_QQ(LSG);
      return Size(Intersection( Kernel(DeterminantOfCharacter(qq.Tchar)) ,
      Kernel(qq.Vchar) ,
      Kernel(qq.sigma) ) ) / Size(LSG!.kerneltau) >1;
      return false;
    else #__HERE__ TODO CHECK
      qq:=SRC_QQ(LSG);
      return Size(Intersection( Kernel(DeterminantOfCharacter(qq.Tchar)) ,
      Kernel(qq.sigma) ) ) >1;
	## return (Size ( Kernel(DeterminantOfCharacter(qq.Tchar))) / Size(LSG!.kerneltau) > 1);
    fi;
  elif LSG!.dim=3
  then
   return false; # not yet... __HERE__ TODO
  else
  return false; ## __HERE__ __TODO__
  fi;
end);


InstallMethod(GroupOrder,"for LSGs", [IsLagSymmetryGroup],
function(LSG)
return Size(SRC_QQ(LSG).group);
end);

InstallMethod(KernelTauOrder,"for LSGs", [IsLagSymmetryGroup],
function(LSG)
return Size( LSG!.kerneltau );
end);
