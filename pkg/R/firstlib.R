.First.lib <- function(lib,pkg)
{
	library.dynam("Survomatic",pkg,lib.loc=lib);
	cat(paste(packageDescription("Survomatic")[c('Package','Version')],collapse=' '));
	cat(" loaded\n");
	data(ctrl,dex,models);
}
