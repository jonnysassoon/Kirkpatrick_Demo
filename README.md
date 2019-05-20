# Kirkpatrick_Demo

This program provides a useful interactive demo of Kirkpatrick's point location algorithm.

### Disclaimer:
For the programmers amongst you that would like to try this demo, and would also like to analyze
my code, I will take this time while I have your attention to point out that this project as
a whole does not exactly reflect best practices in terms of coding style. Hopefully I will
refactor at some point and reduce clutter/redundancies, but such is the nature of final
projects. :)

Here are some general things you need to know for how to use the program

  ## Limitations:
  My implementation is not "complete" in the sense that it only runs correctly
  when subdivisions are created in a certain way. These limitations include:
  * ALL subdivisions must exist inside the initial outer polygon. You cannot provide an initial polygon
    and then re-split the outer plane by clicking on a point along the polygon, and making a new edge
    from that point to somewhere in the outer face. (It should be additionally noted that you cannot
    create an initial polygon, and then add a new polygon - disconnected from the first one - that
    also sits in the outer plane.)
  * There can be no "loose edges" (or vertices of degree 1). Whenever you create a split, you must
    close it in the appropriate way to make a new polygon.
  * You cannot make an entirely new polygon *inside* a pre-existing one. All subsequently created
    polygons must *extend* from a preexisting one, except of course for the intial one (base cases...)
  * When you split, you cannot finish the split by clicking on a vertex that you just created by doing
    *this* split (clicking on another previously split vertex is fine). This includes trying to close
    the split by clicking on the first node of the split.
  All of these limitations should be accounted for in a completed demo of point location as they are 
  legitimate subdivisions.
  
  ## Inputs:
  while there will be a console in the demo to guide you at each step how to input data and progress
  to the next step, there are some general things you should know pertaining to functionality:
  * When you are making the outer polygon, you must hit 'enter' to close it (this is the only time
      you have to do this, see **Assumptions** for why this is).
  * You can hit 'r' at *any time* to completely reset the process and start from scratch.
  ## Assumptions:
  Besides what is listed in the **Limitations** section, there is a key assumption I made for simplicity 
  (not really...I should actually change this at some point):
    When you are making the outer polygon, the program assumes that no points will overlap, i.e. there 
    are no active checks of the sort "have you just clicked on a point you already made?" Consequently,
    you can make the outer polygon look quite awkward by having vertices (which have a radius of 5 pixels)
    that are "overlapping," and you can also likely cause some sort of error if you add a new vertex on
    the same pixel as an already existing vertex. Just avoid this when making the outer polygon.
  ## Bugs:
  While all of the aforementioned might make your usage more annoying, they do not detract from 
  the program running succesfully. But hey, what would a hacky demo be without bugs? Here are a couple of 
  bugs to lookout for:
  * When you're asked to enter a query point, it *must* be inside the outer triangle. For some reason
      I didn't add a check for this (I'm already doing a ton of computation per frame, can't afford to,
      you know...make things work well), and it will cause a nullPointerException and crash the program. 
      Just make sure to only query within the bounds.
  * Sometimes when you split, you could legally finish a split by clicking on a pre-existing vertex, but
      the program continues as though you are continuing to split (as manifested by the program continuing
      to draw the hypothetical edge from the vertex you just clicked on to where your mouse is). I don't know
      why this happens. If you find yourself in this situation, you cannot triangulate the subdivision (even 
      after finishing "that" split) as it will cause a nullPointException. If this happens, you have to start 
      the subdivision process over again by hitting 'r'.
