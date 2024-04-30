
SetPackageInfo( rec(
PackageName := "SymOrb",
Subtitle := "Minimizing orbits",
Version := "0.92",
Date := "15/11/2022",
License := "GPL-3.0",
ArchiveURL :=
          "http://localhost",
ArchiveFormats := ".tgz",
Persons := [
  rec(
    LastName      := "Ferrario",
    FirstNames    := "DL",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "davide.ferrario@unimib.it",
    WWWHome       := "http://www.dlfer.xyz",
    PostalAddress := Concatenation( [
                       "Universit� di Milano-Bicocca\n",
                       "via Cozzi, 52 -- 20125 Milano\n",
                       "Italy" ] ),
    Place         := "Milano",
    Institution   := "University of Milano-Bicocca"
  )
],

Status := "dev",
Dependencies := rec(
  GAP := ">=4.3",
  NeededOtherPackages := [],
  SuggestedOtherPackages := [],
  ExternalConditions := []

),

AvailabilityTest := ReturnTrue,
# AvailabilityTest := function()
#   local path,file;
#     # test for existence of the compiled binary
#     path:=DirectoriesPackagePrograms("example");
#     file:=Filename(path,"hello");
#     if file=fail then
#       Info(InfoWarning,1,
#         "Package ``Example'': The program `hello' is not compiled");
#       Info(InfoWarning,1,
#         "`HelloWorld()' is thus unavailable");
#       Info(InfoWarning,1,
#         "See the installation instructions; ",
#         "type: ?Installing the Example package");
#     fi;
#     # if the hello binary was vital to the package we would return
#     # the following ...
#     #return file<>fail;
#     # since the hello binary is not vital we return ...
#     return true;
#   end,

README_URL :="https://localhost",
PackageInfoURL :="https://localhost",
PackageWWWHome:="https://localhost",
AbstractHTML:="",
PackageWWWHome:="https://localhost",

BannerString := Concatenation(
  "----------------------------------------------------------------\n",
  "Loading  SymOrb", ~.Version, "\n",
  "by ", ~.Persons[1].FirstNames, " ", ~.Persons[1].LastName,
        " (", ~.Persons[1].WWWHome, ")\n",
  "For help, type: ?SymOrb package \n",
  "----------------------------------------------------------------\n" ),

Autoload := false,

##  *Optional*, but recommended: path relative to package root to a file which
##  contains as many tests of the package functionality as sensible.
#TestFile := "tst/testall.g",

##  *Optional*: Here you can list some keyword related to the topic
##  of the package.
# Keywords := ["Smith normal form", "p-adic", "rational matrix inversion"]
Keywords := []

));


