#########################################################################
#########################################################################
####
#### SYMORB
####
#########################################################################
#########################################################################

InstallGlobalFunction(LagSymmetryGroupSRC,
	function (action_type,NOB,dim,kern, rot, ref )
		local LSG;
		LSG:=Objectify(LagSymmetryGroupType,
			rec(
				NOB        := NOB,
				dim	   := dim,
				action_type:= action_type,
				kerneltau  := kern,
				rotation   := rot,
				reflection := ref
			   )
                    );
		if IsValidLSG(LSG) then
			LSG!.action_type:=ActionType(LSG);
			return LSG;
		else
			return fail;
		fi;
	end);

##

InstallGlobalFunction(LagSymmetryGroup,
function (action_type,NOB, kern, rotV, rotS, refV, refS)
## kern subgroup of the direct product...
## should be a group with generators
	local dim;
	dim:=Length(rotV);
	return LagSymmetryGroupSRC(action_type,NOB,dim,kern, Tuple([rotV,rotS]), Tuple([refV,refS]) );
end);


InstallGlobalFunction(LagSymmetryGroupCHARS,
	function (NOB, Tchar, Vchar, sigma )
		local LSG,dim;
		dim:=Degree(Vchar);
		LSG:=Objectify(LagSymmetryGroupCHARSType,
			rec(
				n        := NOB,
				dim	   := dim,
				group := UnderlyingGroup(Tchar),
				Tchar	   := Tchar,
				Vchar	   := Vchar,
				sigma 	:= sigma
			   )
                    );
		return LSG;
	end);

##

##

InstallMethod(MinorbInitElements, "for LagSymmetryGroups", [IsLagSymmetryGroup],
function(LSG)
local result,cyclic_order,i;
if Size( LSG!.kerneltau ) = 1
then
result:= [ [ One (TrivialKerTau(LSG!.dim )) ]  ];
else
result:=[ Elements(LSG!.kerneltau) ];
fi;

cyclic_order:= Size ( Group( Concatenation(GeneratorsOfGroup(LSG!.kerneltau),[LSG!.rotation]) )) / Size(LSG!.kerneltau);
result[2]:=List( [1..cyclic_order], i ->  (LSG!.rotation) ^ i );
if LSG!.action_type > 0
then
result[3]:=LSG!.reflection;
result[4]:=LSG!.rotation * LSG!.reflection;
fi;
return result;
end);






## python float conversion

InstallGlobalFunction(ConvertToString,
function (mat)
    local rows, cols, floatMatrix, i, j;

rows := Length(mat);
    cols := Length(mat[1]);
    floatMatrix := [];
    for i in [1..rows] do
        floatMatrix[i] := [];
        for j in [1..cols] do
            floatMatrix[i][j] := String(mat[i][j]);
        od;
    od;
    return floatMatrix;

end);

##########
##########
##########
##########
InstallGlobalFunction(ConvertToFloatNC,
function (var)
local str,filename,python_script, exec_program, a, tmpstring;
python_script:= "\
import sys,math,cmath,string\n\
\n\
tmpstring = sys.argv[1]\n\
tmpstring = tmpstring.replace(\"^\",\"**\")\n\
tmpstring = tmpstring.replace(\"/\",\"* 1.0 /\")\n\
## sys.stdout.write(tmpstring)\n\
var = eval(tmpstring)+0J \n\
\n\
sys.stdout.write(str(var.real))\n\
";

exec_program := Filename( DirectoriesSystemPrograms(), "python" );
str := ""; a := OutputTextString(str,true);
SetPrintFormattingStatus(a,false);
Process( DirectoryCurrent(), exec_program, InputTextNone(), a, ["-c", python_script, String(var)] );
CloseStream(a);
return str;
end);




############
############
############
############

InstallGlobalFunction(MinorbPrintElement,
function (NOB,dim,tuple,comments)
	local j,k,this_string,output,pperm;
	this_string:=""; output:=OutputTextString(this_string,true);
	SetPrintFormattingStatus(output,false);
	pperm:=Permuted([1..NOB],tuple[2]^(-1));
	for j in [1..NOB] do
        AppendTo(output,String(pperm[j]), " !! ", comments[1], " j=",j,"\n");
    od;
	for j in [1..dim] do
        for k in [1..dim] do
        	AppendTo(output,  ConvertToFloat(tuple[1][j][k])," !! ", comments[2], " (j,k)=(",j,",",k,")\n" );
        od; od;


	CloseStream(output);
	return this_string;
end);

############
############
############
############
InstallGlobalFunction(PrintMatrix,
function(matrix)
local str, a, row, row_number, i, entry ;
str := ""; a := OutputTextString(str,true);
SetPrintFormattingStatus(a,false);
AppendTo(a,"[ ");
for row_number  in  [1..Length(matrix)-1]
do
        row:=matrix[row_number];
        AppendTo(a,"[ ");
        for i in [1..(Length(row)-1)]
                do
                entry:=ConvertToFloatNC(row[i]);
                AppendTo(a, entry, ", ");
                od;
        entry:=ConvertToFloatNC(row[Length(row)]);
        AppendTo(a, entry, " ], ");
od;
row:=matrix[Length(matrix)];
AppendTo(a,"[ ");
for i in [1..(Length(row)-1)]
        do
        entry:=ConvertToFloatNC(row[i]);
        AppendTo(a, entry, ", ");
        od;
entry:=ConvertToFloatNC(row[Length(row)]);

AppendTo(a, entry, " ] ]");
CloseStream(a);
return str;
end);

###########
###########
###########
###########

InstallGlobalFunction(MinorbDisplayElement,
function (NOB,dim,tuple,comments)
	local j,k,this_string,output,pperm;
	this_string:=""; output:=OutputTextString(this_string,true);
	SetPrintFormattingStatus(output,false);
	pperm:=Permuted([1..NOB],tuple[2]^(-1));


	for j in [1..NOB]
                do
		AppendTo(output,
"x_", j, "(",comments[1],"t",comments[2],") = " ,
PrintMatrix(tuple[1]), " * x_", pperm[j],
"(t)\n");
                 od;

	CloseStream(output);
	return this_string;
end);


########
########
########
########

InstallMethod(MinorbInitString, "for LagSymmetryGroups", [IsLagSymmetryGroup],
function(LSG)
local tmplist,init_string,output,i,elm,comments,cyclic_order;
init_string:="";
output:=OutputTextString(init_string,true);
SetPrintFormattingStatus(output,false);
tmplist:=MinorbInitElements(LSG);

AppendTo(output,"\
!!!======================================================================\n\
!!! generated configuration file\n\
!!!======================================================================\n");

# &CONFIG .. /\n later
## omegax, omegay,omegaz later

##computed first... LSG!.action_type:=ActionType(LSG);
AppendTo(output,String(LSG!.action_type), " !! action_type\n" );
AppendTo(output,String(LSG!.NOB), " !! NOB\n" );
AppendTo(output,String(LSG!.dim), " !! dim\n" );

# for i in [1..LSG!.NOB]
# do
# 	AppendTo(output,String(1), "  !! m(", String(i), ")\n");
# od;

AppendTo(output, String(Size(LSG!.kerneltau)), "  !! sizekerT\n");

comments:=[ "perm", "matrix" ];

for elm in tmplist[1]
do
	AppendTo(output,MinorbPrintElement(LSG!.NOB,LSG!.dim,elm,comments));
od;

cyclic_order:=Length(tmplist[2]);
AppendTo(output, cyclic_order, "  !! cyclic_order\n");

comments:= [ "cyclic_perm", "cyclic_matrix" ];
for elm in tmplist[2]
do
	AppendTo(output,MinorbPrintElement(LSG!.NOB,LSG!.dim,elm,comments));
od;

if LSG!.action_type>0
then
    comments := [ "first_perm", "first_matrix" ];
    elm:=tmplist[3];
    AppendTo(output,MinorbPrintElement(LSG!.NOB,LSG!.dim,elm,comments));
    comments := [ "second_perm", "second_matrix" ];
    elm:=tmplist[4];
    AppendTo(output,MinorbPrintElement(LSG!.NOB,LSG!.dim,elm,comments));
fi;

CloseStream(output);
return init_string;
end);
########
########
########
########

InstallMethod(MinorbInfoString, "for LagSymmetryGroups", [IsLagSymmetryGroup],
function(LSG)
local tmplist,info_string,output,i,elm,comments,cyclic_order,small_gens,
path,date;
info_string:="";
output:=OutputTextString(info_string,true);
SetPrintFormattingStatus(output,false);

AppendTo(output,"% SYMORB version : ", SYMORB_VERSION, "\n");
path:=DirectoriesSystemPrograms();
date:=Filename(path,"date");
AppendTo(output,"% date        : ");
Process(DirectoryCurrent(), date, InputTextUser(), output , []);
AppendTo(output,"% on          : ");
date:=Filename(path,"uname");
#__ Process(DirectoryCurrent(), date, InputTextUser(), output , ["--machine", "--sysname", "--processor"]);
Process(DirectoryCurrent(), date, InputTextUser(), output , ["-m", "-s", "-p"]);


##
# begin_attributes
##
AppendTo(output,"@ GroupOrder: ", String(GroupOrder(LSG)),"\n");
AppendTo(output,"@ KernelTauOrder: ", String(KernelTauOrder(LSG)),"\n");
AppendTo(output,"@ ActionType: ", String(ActionType(LSG)),"\n");
## AppendTo(output,"@ IsTypeR: ", String(IsTypeR(LSG)),"\n");
AppendTo(output,"@ IsCoercive: ", String(IsCoercive(LSG)),"\n");
## AppendTo(output,"@ IsRedundant: ", String(IsRedundant(LSG)),"\n");
## AppendTo(output,"@ HasAlwaysCollisions: ", String(HasAlwaysCollisions(LSG)),"\n");
## AppendTo(output,"@ TransitiveDecomposition: ", String(TransitiveDecomposition(LSG)),"\n");
## AppendTo(output,"@ TypeRDirections: ", String(TypeRDirections(LSG)),"\n");

##
# end_attributes
##
cyclic_order:= Size ( Group( Concatenation(GeneratorsOfGroup(LSG!.kerneltau),[LSG!.rotation]) )) / Size(LSG!.kerneltau);

if Size(LSG!.kerneltau)>1 then
small_gens:=GeneratorsOfGroup(LSG!.kerneltau);

for elm in small_gens
do
comments:=[ "" , "" ];
AppendTo(output,MinorbDisplayElement(LSG!.NOB,LSG!.dim, elm ,comments  ));
od;
fi;

## Rotation
if cyclic_order>1 then
    comments:=[ "", Concatenation("+T/",String(cyclic_order))  ];
    AppendTo(output,MinorbDisplayElement(LSG!.NOB,LSG!.dim, LSG!.rotation ,comments  ));
fi;

## Reflection
if LSG!.action_type>0
then
comments:=[ "-", ""  ];
AppendTo(output,MinorbDisplayElement(LSG!.NOB,LSG!.dim, LSG!.reflection ,comments  ));
fi;


CloseStream(output);
return info_string;
end);


#########
#########
#########
#########



InstallGlobalFunction(MakeMinorbSymFile,
function (filename,LSG)
local symfile,infofile,fdsymfile,fdinfofile;
symfile:=Concatenation(filename,".sym");
infofile:=Concatenation(filename,".info");
fdsymfile:=OutputTextFile(symfile,false);
fdinfofile:=OutputTextFile(infofile,false);
SetPrintFormattingStatus(fdsymfile,false);
SetPrintFormattingStatus(fdinfofile,false);
PrintTo(fdsymfile,MinorbInitString(LSG));
CloseStream(fdsymfile);
Print("file ", symfile , " created!\n");
PrintTo(fdinfofile,MinorbInfoString(LSG));
CloseStream(fdinfofile);
Print("file ", infofile , " created!\n");
return 0;
end);

##
##
##
InstallGlobalFunction(TrivialKerTau,
function(dim)
local trivialmat,i;
trivialmat:=DiagonalMat(List([1..dim], i -> 1));
return GroupWithGenerators( [ Tuple( [trivialmat,()]) ] );
end);


InstallGlobalFunction(CyclicLagSymmetryGroupTK,
function (NOB,rotV, rotS, refV, refS)
## kern subgroup of the direct product...
## should be a group with generators
	local dim,action_type,kern;
	dim:=Length(rotV);
	action_type:=0;
	kern:=TrivialKerTau(dim);
	return LagSymmetryGroupSRC(action_type,NOB,dim,kern, Tuple([rotV,rotS]), Tuple([refV,refS]) );
end);


InstallGlobalFunction(LagSymmetryGroupTK,
function (action_type,NOB,rotV, rotS, refV, refS)
## kern subgroup of the direct product...
## should be a group with generators
	local dim,kern;
	dim:=Length(rotV);
	kern:=TrivialKerTau(dim);
	return LagSymmetryGroupSRC(action_type,NOB,dim,kern, Tuple([rotV,rotS]), Tuple([refV,refS]) );
end);


InstallGlobalFunction(SRCProduct,
function(action_type,list)
local l,i;
l:=Length(list);
if action_type=0 then
	return DiagonalMat(List([1..l], i -> list[i]   ));
else
	return [
		DiagonalMat(List([1..l], i -> list[i][1]   )),
		DiagonalMat(List([1..l], i -> list[i][2]   ))
	       ];
fi;

end);

InstallGlobalFunction(AllRhosTK,
function(action_type,dim)
	local mylist,result;
	if (action_type = 0) then
		mylist:=[  1   ,   -1     ];
	else
		mylist:= [ [ 1,1 ],
			   [ 1, -1 ],
			   [ -1,1 ],
			   [ -1,-1]
			];
	fi;

	if (dim=0)
	then
		Error("dimension wrong\n");
	else
		return List( UnorderedTuples( mylist, dim ), x ->
			SRCProduct(action_type,x)
			);
	fi;
end);

InstallGlobalFunction(SRCSum,
function(action_type,list)
local l,i;
l:=Length(list);
if action_type=0 then
	return 0;
else
	return [0 , 0 ];
fi;

end);


InstallGlobalFunction(AllTransitiveSigmasTK,
function(action_type,NOB)
local i;
if action_type=0 then
	return [ PermList( List( [1..NOB] , i ->  i mod NOB +1  ))] ;
else
    if NOB=1 then
	return [ [(),()] ];
    elif NOB=2 then
        return [ [(),(1,2)],
                 [(1,2),(1,2)],
		 [(1,2),() ]
		];
    elif NOB=3 then
	return [ [ (1,2,3), (1,2) ] ];
    else
	Error( "not yet implemented\n");
	return 0;
    fi;

fi;
end);

InstallGlobalFunction(AllSigmasTK,
function(action_type,NOB)
local p,res,tmpdata,i;
tmpdata:=List([1..NOB], i -> AllTransitiveSigmasTK(action_type,i));
res:=[];
if NOB <> 3 then
    Error("not yet implemented for NOB <> 3\n");
    return;
fi;
## 1+1+1
Append(res,tmpdata[1]);
## 2+1
Append(res,tmpdata[2]);
## 3
Append(res,tmpdata[3]);
return res;
end);

InstallGlobalFunction(AllLSGTK,
function(action_type,NOB,dim)
local allsigmas, allrhos,rho,sigma, res, rotV, rotS, refV,refS;
allsigmas:=AllSigmasTK(action_type,NOB);
allrhos:=AllRhosTK(action_type,dim);
res:=[];
for sigma in allsigmas do
	for rho in allrhos do
	if action_type=0 then
	    rotV:=rho;
	    rotS:=sigma;
	    refV:=-1;
	    refS:=-1;
	else
	    rotV:=rho[1];
	    rotS:=sigma[1];
	    refV:=rho[2];
	    refS:=sigma[2];
	fi;
	Append(res, [
	LagSymmetryGroupTK(action_type,NOB,rotV,rotS,refV,refS)
	 ]);
	od;
od;
return res;
end);



## todo: classification AllLagSymmetryGroups(Size,nnn ... )
## todo: IsTransitive
## todo: IsTrivial
## todo: IsAlwaysCollision
## todo: ....

## todo: kernel building (GroupWithGenerators)

