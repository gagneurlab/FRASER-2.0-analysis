
theme_manuscript <- function(fig_font="Arial", fig_font_size=8){ 
    font <- fig_font   #assign font family up front
    
    theme_pubr() %+replace%    #replace elements we want to change
        
        theme(
            
            #grid elements
            panel.grid.major = element_blank(),    #strip major gridlines
            panel.grid.minor = element_blank(),    #strip minor gridlines
            # axis.ticks = element_blank(),          #strip axis ticks
            
            #since theme_minimal() already strips axis lines, 
            #we don't need to do that again
            
            #text elements
            text = element_text(          
                family = font,            
                size = fig_font_size),
            
            # axis labels
            axis.title=element_text(face="bold"), 
            
            # legend
            legend.position="top"
        )
}
