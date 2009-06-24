.First.lib <- function(lib,pkg)
{
	library.dynam("Survomatic",pkg,lib.loc=lib);
	cat("Survomatic 0.1.2.4 loaded\n");
	data(ctrl,dex,models);
}
