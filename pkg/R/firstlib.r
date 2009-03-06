.First.lib <- function(lib,pkg)
{
	library.dynam("Survomatic",pkg,lib.loc=lib);
	cat("Survomatic 0.1-1 loaded\n");
	data(ctrl,dex,exgm,exg,exlm,exl);
}
