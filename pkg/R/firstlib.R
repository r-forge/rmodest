.First.lib <- function(lib,pkg)
{
	library.dynam("Survomatic",pkg,lib.loc=lib);
	cat(paste(packageDescription("Survomatic")[c('Package','Version')],collapse=' '));
	cat(" loaded\n");
	data(ctrl,dex,modelinfo,modelpars,envir=environment());
	options(ctrl=ctrl);
	options(dex=dex);
	options(dmodels=dmodels);
	options(modelpars=modelpars);
}