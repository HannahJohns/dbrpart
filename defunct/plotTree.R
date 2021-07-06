#' Combine template and node plots to create something
#' approximating a conventional CART tree
#'
#'
#'
plotTree <- function(template,plots,options = NULL)
{
  if(is.null(options))
  {
    options <- list(tileWidth = 1000,
                    tileHeight=400,
                    tileSize = 0.1,
                    plotTerminal = TRUE)
  }

  # Get width/height (% of canvas) for each tile from inputs
  width <- options$tileSize
  height <- options$tileHeight/options$tileWidth*options$tileSize

  # This should eventually be removed
  if(options$plotTerminal)
  {
    warning("Terminal Plotting is not available yet. Setting to FALSE")
    options$plotTerminal <- FALSE
  }

  # If we're not plotting terminal nodes, template's positioning
  # will be incorrect. Recalculate shift/scale for template
  if(!options$plotTerminal)
  {
    xmin <- min(template[which(template$terminal==0),"x"])
    ymin <- min(template[which(template$terminal==0),"y"])

    template$x <- template$x-xmin
    template$xInEnd <- template$xInEnd-xmin
    template$xOutEnd <- template$xOutEnd-xmin

    template$y <- template$y-ymin
    template$yEnd <- template$yEnd-ymin

    # Scale x and y separately so top right of tile is in

    xmax <- max(template[which(template$terminal==0),"x"])+width

    template$x <- template$x/xmax
    template$xInEnd <- template$xInEnd/xmax
    template$xOutEnd <- template$xOutEnd/xmax

    ymax <- max(template[which(template$terminal==0),"y"])+height

    template$y <- template$y/ymax
    template$yEnd <- template$yEnd/ymax


  }




  # Template is duplicated for each class,
  # which isn't needed here.
  # Drop duplicates
  template <- unique(template[,-(2:3)])

  # Set up the canvas
  ggp <- cowplot::ggdraw()

  for (node in names(plots))
  {
    # What happens at this node
    # depends on if we said we wanted to plot
    # terminal nodes or not

    if(options$plotTerminal | template[which(as.character(template$label)==node),"terminal"]==0)
    {
      # If this node isn't terminal, or if we said we should plot terminal nodes
      # Prepare the tile for this node

      print(node)
      tile <- magick::image_graph(width=options$tileWidth,height=options$tileHeight)
      dev.set(which = dev.list()["magick"])
        plot.nodePlot(plots[[node]])
      dev.off(which = dev.list()["magick"])

    }

    # Drawing the stems connecting nodes is a bit more complicated
    # and requires us to look forward down the tree.
    # First, check that this node isn't terminal
    if(template[which(as.character(template$label)==node),"terminal"]==0)
    {

      inNodeTerminal <- template[which(template$label==
                                      template[which(as.character(template$label)==node),
                                      "inLabel"]
                                      ),
                                      "terminal"]


      if(options$plotTerminal |
         unique(ifelse(length(inNodeTerminal)>0,inNodeTerminal,1))==0
      )
      {
        ggp <- ggp + cowplot::draw_line(x=unlist(template[which(as.character(template$label)==node),c("x","xInEnd")]) + width*0.5,
                                        y=unlist(template[which(as.character(template$label)==node),c("y","yEnd")]) + height*0.5,
                                        color="black"
        )
      }



      outNodeTerminal <- template[which(template$label==                                      # The node corresponding
                                        template[which(as.character(template$label)==node),   # to the in group for this split
                                        "outLabel"]
                                        ),                                          # is not terminal
                                  "terminal"]
      if(options$plotTerminal |                                               # If we plot terminal nodes, or
         unique(ifelse(length(outNodeTerminal)>0,outNodeTerminal,1))==0
      )
      {

        ggp <- ggp + cowplot::draw_line(x=unlist(template[which(as.character(template$label)==node),c("x","xOutEnd")]) + width*0.5,
                               y=unlist(template[which(as.character(template$label)==node),c("y","yEnd")]) + height*0.5,
                               color="black"
        )

      }
    }



    if(options$plotTerminal | template[which(as.character(template$label)==node),"terminal"]==0)
    {
      # If this node isn't terminal, or if we said we should plot terminal nodes
      # Add the tile to the canvas

      ggp <- ggp + cowplot::draw_image(tile,
                                       width=width,
                                       height=height,
                                       x=template[which(as.character(template$label)==node),"x"],
                                       y=template[which(as.character(template$label)==node),"y"])

    }


  }


  return(ggp)
}
