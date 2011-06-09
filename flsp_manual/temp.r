
fitted <- as.data.frame(nzrl@fitted_index)
fitted$qname = "fitted_index"
index <- rbind(as.data.frame(nzrl@index),fitted)
xyplot(data ~ year, group=qname, data=index, type="b",auto.key=TRUE)

# If multiple indices
nz <- nzrl
nz@index[[2]] <- nz@index[[1]]
names(nz@index)[2] <- "index2"
nz@fitted_index[[2]] <- nz@fitted_index[[1]]
names(nz@fitted_index)[2] <- "index2"
nz@residuals_index[[2]] <- nz@residuals_index[[1]]
names(nz@residuals_index)[2] <- "index2"

# Plot the index
fitted <- cbind(as.data.frame(nz@fitted_index),type="fitted")
index <- cbind(as.data.frame(nz@index),type="index")
index <- rbind(index,fitted)
xyplot(data ~ year | qname, group=type, data=index, type="b",auto.key=TRUE)

# Plot the residuals
residuals <- as.data.frame(nz@residuals_index)
#xyplot(data ~ year | qname, data=residuals)
xyplot(data ~ year | qname, data=residuals,panel=function(x,y)
	{panel.xyplot(x,y)
	panel.loess(x,y,span=1) })


# Table of results




