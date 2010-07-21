.First.lib <- function(lib,pkg)
{
	library.dynam("Survomatic",pkg,lib.loc=lib);
	cat(paste(packageDescription("Survomatic")[c('Package','Version')],collapse=' '));
	cat(" loaded\n");
	data(ctrl,dex,modeldeps,modelpars,envir=environment());
	options(ctrl=ctrl);
	options(dex=dex);
	options(modeldeps=modeldeps);
	options(modelpars=modelpars);
}