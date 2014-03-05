.onAttach <- function(libname, pkgname) {

  banner.text <- paste("\nWelcome to the glmnetGLR package, version ",
                       packageDescription("glmnetGLR", fields="Version"), ".\n\n",
                       "Commented source code can be found in\n",
                       path.package(package="glmnetGLR"), "/SourceCode.\n", sep="")
  
  citation.text <- "\nPlease cite the following reference: 
\nAmidan BG, Orton DJ, LaMarche BL, Monroe ME, Moore RJ, 
Venzin AM, Smith RD, Sego LH, Payne SH, Tardiff MF. 
Signatures for Mass Spectrometry Data Quality.
(submitted for publication)"

  packageStartupMessage(banner.text)
  packageStartupMessage(citation.text)
  
}


