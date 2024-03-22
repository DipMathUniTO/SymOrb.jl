#########################################################################
#########################################################################
####
#### SYMORB
####
#########################################################################
#########################################################################

# KRH decomposition: sum, product, ...
## UnbindGlobal("LagSymmetryGroupFam");
LagSymmetryGroupFamily:=NewFamily("LagSymmetryGroupFam");

## UnbindGlobal( "IsLagSymmetryGroup" );
DeclareRepresentation("IsLagSymmetryGroup",IsComponentObjectRep, 
[
"NOB","dim","action_type",
"kerneltau",
"rotation",
"reflection"
]
); 



## UnbindGlobal( "IsLagSymmetryGroupCHARS" );
DeclareRepresentation("IsLagSymmetryGroupCHARS",IsComponentObjectRep, 
[ "n","dim", "group", "Tchar", "Vchar", "sigma" ]
); 
## UnbindGlobal("LagSymmetryGroupCHARSType");
LagSymmetryGroupCHARSType:=NewType(LagSymmetryGroupFamily, IsLagSymmetryGroupCHARS);

InstallMethod(PrintObj,"for LSGCHARSs",
  [IsLagSymmetryGroupCHARS],
function(LSG)
  Print("LagSymmetryGroupCHARS(NOB=",
	LSG!.n,
	", dim=", LSG!.dim,
	", Tchar=", LSG!.Tchar,
	", Vchar=", LSG!.Vchar,
	", sigma=", LSG!.sigma,
	"\n"
 );
end);

DeclareGlobalFunction("LagSymmetryGroupCHARS");








InstallMethod(PrintObj,"for LSGs",
  [IsLagSymmetryGroup],
function(LSG)
 if LSG!.action_type<>0 then
  Print("LagSymmetryGroup(NOB=",
	LSG!.NOB,
	", dim=",LSG!.dim,
	", action_type=", LSG!.action_type,
	", rotation=", LSG!.rotation,
	", reflection=", LSG!.reflection,
	"\n"
 );
else
  Print("LagSymmetryGroup(NOB=",
	LSG!.NOB,
	", dim=",LSG!.dim,
	", action_type=", LSG!.action_type,
	", rotation=", LSG!.rotation,
	"\n"
 );
fi;

end);



## UnbindGlobal("LagSymmetryGroupType");
LagSymmetryGroupType:=NewType(LagSymmetryGroupFamily, IsLagSymmetryGroup);


## UnbindGlobal("LagSymmetryGroupSRC");
DeclareGlobalFunction("LagSymmetryGroupSRC");


##

## UnbindGlobal("LagSymmetryGroup");
DeclareGlobalFunction("LagSymmetryGroup");


##
## UnbindGlobal("MinorbInitElements");
DeclareAttribute("MinorbInitElements",IsLagSymmetryGroup);

##
## UnbindGlobal("ConvertToFloat");
DeclareGlobalFunction("ConvertToFloat");

## UnbindGlobal("ConvertToFloatNC");
DeclareGlobalFunction("ConvertToFloatNC");



##
## UnbindGlobal("MinorbPrintElement");
DeclareGlobalFunction("MinorbPrintElement");

## UnbindGlobal("PrintMatrix");
DeclareGlobalFunction("PrintMatrix");

## UnbindGlobal("MinorbDisplayElement");
DeclareGlobalFunction("MinorbDisplayElement");

##
## UnbindGlobal("MinorbInitString");
DeclareAttribute("MinorbInitString",IsLagSymmetryGroup);

## UnbindGlobal("MinorbInfoString");
DeclareAttribute("MinorbInfoString",IsLagSymmetryGroup);

##
## UnbindGlobal("MakeMinorbSymFile");
DeclareGlobalFunction("MakeMinorbSymFile");

DeclareGlobalFunction("TrivialKerTau");

## UnbindGlobal("CyclicLagSymmetryGroupTK");
DeclareGlobalFunction("CyclicLagSymmetryGroupTK");

## UnbindGlobal("LagSymmetryGroupTK");
DeclareGlobalFunction("LagSymmetryGroupTK");

## UnbindGlobal("SRCProduct");
DeclareGlobalFunction("SRCProduct");

## UnbindGlobal("SRCProduct");
DeclareGlobalFunction("SRCSum");

## UnbindGlobal("AllRhosTK");
DeclareGlobalFunction("AllRhosTK");

## UnbindGlobal("AllTransitiveSigmasTK");
DeclareGlobalFunction("AllTransitiveSigmasTK");

## UnbindGlobal("AllSigmasTK");
DeclareGlobalFunction("AllSigmasTK");

## UnbindGlobal("AllLSGTK");
DeclareGlobalFunction("AllLSGTK");

