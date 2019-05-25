import java.util.*;
import java.lang.Math;

int TREEXSTART = 610;
int CONSHEIGHT = 190;
int WORLDWIDTH = 1000;
int WORLDHEIGHT = 750;
Point LEFT = new Point(20, 720, -3);
Point RIGHT = new Point(600, 720, -2);
Point TOP = new Point(310, 200, -1); // vertices of the triangle
ArrayList<Point> outerPolygon, currentPolygon;
ArrayList<ArrayList<Point>> completedPolygons;
HashSet<MapTriangle> triangulation;
HashSet<Point> indSet;
MapTriangle TRI; // parent triangle, always keep ref to this
MapTriangle CURRTRI; // triangle we're in
MapTriangle NEXTTRI; // the next triangle in the tree that the query point is inside
boolean completedFirstPolygon, isSplitting, isClockwise;
int polygonToSplit, idNum, result, triNum; // when we split, we need to know the index of the polygon we're splitting
Point QueryPoint;
HashSet<MapTriangle> allTris; // have a running set of ALL triangles for printing
HashMap<Integer, Integer> rankAtDepth;
int treeDepth;

// TODO: if we query outside the triangle, it breaks.
// TODO: add labels to polygons
// TODO: make the tree a bit better
  // TODO: add tree level next to each level
  // TODO: figure out appropriate font of text

void setup() {
  size(1000, 750);
  completedFirstPolygon = isSplitting = isClockwise = false;
  outerPolygon = new ArrayList<Point>();
  currentPolygon = new ArrayList<Point>();
  indSet = new HashSet<Point>();
  completedPolygons = new ArrayList<ArrayList<Point>>();
  triangulation = new HashSet<MapTriangle>();
  polygonToSplit = -1; // map -1 to the outer face (between triangle and polygon)
  idNum = triNum = 0; // unique identifier for each node and triangle
  result = -2;
  QueryPoint = null;
  TRI = CURRTRI = NEXTTRI = null;
  allTris = new HashSet<MapTriangle>();
  rankAtDepth = new HashMap<Integer, Integer>();
  treeDepth = 0;
}

void draw() {
  background(0, 0, 255);
  fill(0);
  if (result != -1) fill(128, 128, 255);
  else fill(0, 255, 0);
  triangle(LEFT.x, LEFT.y, TOP.x, TOP.y, RIGHT.x, RIGHT.y);
  fill(0);
  ellipse(LEFT.x, LEFT.y, 15, 15); 
  ellipse(TOP.x, TOP.y, 15, 15);
  ellipse(RIGHT.x, RIGHT.y, 15, 15);
  Point mouseP = new Point(mouseX, mouseY);
  drawConsole();
  if (triangulation.size() > 1 || TRI !=null) {
    drawTriTree();
  }
  if (result < -1 && (triangulation.size() > 1 || CURRTRI != null)) { // work on this
    drawTriangulation();
    if (QueryPoint != null) {
      fill(0, 255, 0);
      ellipse(QueryPoint.x, QueryPoint.y, 5, 5);
      fill(0);
      for (int i = 0; i < completedPolygons.size(); i++) {
        drawPolygon(completedPolygons.get(i), true, i);
      }
      if (CURRTRI.children.size() == 0) {
        fill(255, 255, 0);
        beginShape();
        vertex(CURRTRI.verts[0].x, CURRTRI.verts[0].y);
        vertex(CURRTRI.verts[1].x, CURRTRI.verts[1].y);
        vertex(CURRTRI.verts[2].x, CURRTRI.verts[2].y);
        endShape(CLOSE);
      } else{
        strokeWeight(3);
        stroke(255, 0, 0);
        line(NEXTTRI.verts[0].x, NEXTTRI.verts[0].y, NEXTTRI.verts[1].x, NEXTTRI.verts[1].y);
        line(NEXTTRI.verts[1].x, NEXTTRI.verts[1].y, NEXTTRI.verts[2].x, NEXTTRI.verts[2].y);
        line(NEXTTRI.verts[2].x, NEXTTRI.verts[2].y, NEXTTRI.verts[0].x, NEXTTRI.verts[0].y);
        strokeWeight(1);
        stroke(0);
      }
    }
  }
  else {
    if (result == -1) {
      fill(128, 128, 255);
      beginShape();
      for (int i = 0; i < outerPolygon.size(); i++) {
        vertex(outerPolygon.get(i).x, outerPolygon.get(i).y);
      }
      endShape(CLOSE);
    }
    if (!completedFirstPolygon && outerPolygon.size() > 0) {
      drawPolygon(outerPolygon, false, -1);
    } 
    else{ // draw all polygons
      for (int i = 0; i < completedPolygons.size(); i++) {
        if (i == result) {
          fill(0, 255, 0);
          beginShape();
          for (Point p : completedPolygons.get(i)) {
            vertex(p.x, p.y);
          }
          endShape(CLOSE);
        }
        drawPolygon(completedPolygons.get(i), true, i);
      }
      if (isSplitting) {
        //Point prev = currentPolygon.get(currentPolygon.size() - 1);
        //if (noIntersects(mouseP, currentPolygon, completedPolygons)) stroke(0);
        //else stroke(255, 0, 0);
        //line(mouseP.x, mouseP.y, prev.x, prev.y);
      }
      stroke(0);
      drawPolygon(currentPolygon, false, -1);
    }
    if (TRI != null && result < -1) {
      fill(0, 255, 0);
      ellipse(mouseP.x, mouseP.y, 5, 5);
    }
  }
  if (!completedFirstPolygon) { // draw outerPolygon
    if (outerPolygon.size() == 0) writeToConsole("Click in the triangle to start drawing the polygon");
    else if (outerPolygon.size() == 1) writeToConsole("Drawing the outer polygon...");
    else writeToConsole("Drawing the outer polygon...hit 'enter' to finish");
  } 
  else if (!isSplitting && triangulation.size() == 0 && TRI == null) {
    writeToConsole("Click on a vertex if you want to split, or hit 'tab' to triangulate");
  } 
  else if (isSplitting) {
    writeToConsole("Now splitting...Click on a vertex (not on your split) to finish");
  }
  else if (TRI == null) {
    if (indSet.size() == 0) writeToConsole("Building Graph: Hit 'backspace' to get independent set");
    else writeToConsole("Building Graph: Hit 'backspace' to delete the next point");
  } 
  else if (QueryPoint == null) {
    writeToConsole("Graph has been built. Click to enter query point");
  } 
  else if (result  < -1) {
    if (triangulation.size() < 1) writeToConsole("Reached the bottom of the graph...Click anywhere to continue");
    else writeToConsole("Now traversing the graph...Click anywhere to continue");
  } 
  else {
    writeToConsole("Found polygon " + result + ". Hit 'n' to enter new query point, 'r' to start over"); 
  }
}

void mouseClicked() { // TODO: make this more concise
  Point point = new Point(mouseX, mouseY, idNum);
  if (TRI != null) {
    if (QueryPoint == null) QueryPoint = point;
    if (CURRTRI == null) { // when do you reset this to null? when you finish querying
      CURRTRI = TRI; // current is the top
      triangulation = CURRTRI.children;
      NEXTTRI = CURRTRI.getNextTri(QueryPoint);
    } 
    else {
      if (triangulation.size() == 0) {
        result = CURRTRI.faceNum;
      }
      else {
        CURRTRI = NEXTTRI;
        triangulation = CURRTRI.children;
        NEXTTRI = CURRTRI.getNextTri(QueryPoint);
      }
    }
  }
  if (!completedFirstPolygon) { // Note: the "outer" polygon will change after a split
    if (outerPolygon.size() > 0) {
      Point prev = outerPolygon.get(outerPolygon.size() - 1);
      if (!isInsideTriangle(point, LEFT, TOP, RIGHT) || !noIntersects(point, outerPolygon, completedPolygons)) {
        stroke(255, 0, 0);
        line(point.x, point.y, prev.x, prev.y);
        fill(0);
        stroke(0);
      }
    }
    addToOuterPolygon(point); // Note: this assumes that the new point will *never* be the same as a previous point  
  } if (!isSplitting) { // check if it's a valid split
    Point splitPoint = getSplitSpot(point);
    if (splitPoint != null) {
      addPointToPolygon(splitPoint, currentPolygon);
      isSplitting = true;
    }
  } 
  else {
    Point otherSplitPoint = getSplitSpot(point);
    if (otherSplitPoint == null) otherSplitPoint = point;
    if ( !noIntersects(otherSplitPoint, currentPolygon, completedPolygons) ) {
      Point prev = currentPolygon.get(currentPolygon.size() - 1);
      stroke(255, 0, 0);
      line(point.x, point.y, prev.x, prev.y);
      stroke(0);
      return;
    }
    if (polygonToSplit == -1) { // we have not yet defined the polygons we're splitting
      if (currentPolygon.size() == 1 && !otherSplitPoint.equals(point)) { // a direct split
        // we can't define that point to be on a unique polygon
        // so we can take the midway point
        Point prev = currentPolygon.get(0);
        Point dummyP = new Point( (prev.x + otherSplitPoint.x)/2, (prev.y + otherSplitPoint.y)/2);
        polygonToSplit = insidePolygon(dummyP, completedPolygons);
      } else {
        polygonToSplit = insidePolygon(otherSplitPoint, completedPolygons);
      }
      if (polygonToSplit < 0) return; // we weren't inside one of the polygons
    }
    currentPolygon.add(otherSplitPoint);
    if (!otherSplitPoint.equals(point)) {
      // we snapped to some other vertex that we know does not cause an intersect and 
      // doesn't equal a vertex that's already in the split (see segmentsIntersect logic)
      // i.e. we're finishing a split
      makeSplit(currentPolygon, polygonToSplit);
      currentPolygon = new ArrayList<Point>(); // reset
      isSplitting = false;
      polygonToSplit = -1;
    } else {
      drawPoint(otherSplitPoint);
      idNum++;
    }
  }
}

void keyPressed() { // TODO: have trigger for triangulation
  if (key == 'r') {
    setup();
    return;
  }
  if (key == 'n') {
    QueryPoint = null;
    result = -2;
    NEXTTRI = CURRTRI = null;
    return;
  }
  if (key == ENTER && completedFirstPolygon == false) {
    if ( // illegal end point 
      outerPolygon.size() < 3 ||
      segmentIntersectsWithPolygon(
        outerPolygon.get(outerPolygon.size()-1), 
        outerPolygon.get(0), 
        outerPolygon, false, true)
    ) return;
    isClockwise = orientationIsClockwise(outerPolygon);
    completedPolygons.add(outerPolygon);
    completedFirstPolygon = true;
  }
  if (key == TAB && completedFirstPolygon == true && currentPolygon.size() == 0) {
    triangulate(completedPolygons);
  }
  if (key == BACKSPACE && triangulation.size() > 1) {
    HashMap<Point, HashMap<String, Object>> lookupTable;
    lookupTable= new HashMap<Point, HashMap<String, Object>>();
    makeTable(lookupTable, triangulation); // this is wildly inefficient
    HashSet<Point> marked = new HashSet<Point>();
    marked.add(LEFT);
    marked.add(RIGHT);
    marked.add(TOP);
    if (indSet.size() == 0) getIndSet(lookupTable, marked);
    else {
      removeIndPoint(lookupTable);
    }
  }
  if (triangulation.size() == 1) {
    Iterator iter = triangulation.iterator();
    MapTriangle first = (MapTriangle) iter.next();
    TRI = first;
  }
}

int getRankAtDepth(MapTriangle mt) {
  int rank = 0;
  if (rankAtDepth.containsKey(mt.depth)) {
    rank = rankAtDepth.get(mt.depth)+1;
    rankAtDepth.put(mt.depth, rank);
  } else { // we have a new depth
    rankAtDepth.put(mt.depth, rank);
    treeDepth = mt.depth;
  }
  return rank;
}

void putTriInTree(MapTriangle mt) {
  if (allTris.contains(mt)) {
    return;
  }
  mt.lvlRank = getRankAtDepth(mt);
  allTris.add(mt);
}

void drawTri(Triangle tri) {
  if (indSet.contains(tri.verts[0])) {
    fill(255, 0, 0);
    ellipse(tri.verts[0].x, tri.verts[0].y, 10, 10);
  }
  else drawPoint(tri.verts[0]);
  if (indSet.contains(tri.verts[1])) {
    fill(255, 0, 0);
    ellipse(tri.verts[1].x, tri.verts[1].y, 10, 10);
  }
  else drawPoint(tri.verts[1]);
  if (indSet.contains(tri.verts[2])) {
    fill(255, 0, 0);
    ellipse(tri.verts[2].x, tri.verts[2].y, 10, 10);
  }
  else drawPoint(tri.verts[2]);
  drawLineSegment(tri.verts[0], tri.verts[1], 128);
  drawLineSegment(tri.verts[1], tri.verts[2], 128);
  drawLineSegment(tri.verts[2], tri.verts[0], 128);
}

void drawTriangulation() {
  for (MapTriangle tri : triangulation) {
    drawTri(tri);
  }
}

void drawPolygon(ArrayList<Point> pgon, boolean finished, int label) { //TODO: draw label inside polygon
  for (int j = 0; j < pgon.size(); j++) {
    drawPoint(pgon.get(j));
    if ( j < pgon.size()-1 ) drawLineSegment(pgon.get(j), pgon.get(j+1));
    else if (finished) drawLineSegment(pgon.get(0), pgon.get(pgon.size()-1));
  }
  if (!finished) return;
  Point minX = pgon.get(getMinXindex(pgon));
  Point maxX = pgon.get(getMaxXindex(pgon));
  Point minY = pgon.get(getMinYindex(pgon));
  Point maxY = pgon.get(getMaxYindex(pgon));
  text(label, (minX.x + maxX.x)/2, (minY.y + maxY.y)/2);
}

float getTextWidth(int d) {
  int trisAtDepth = rankAtDepth.get(d);
  float textWidth = (1000-TREEXSTART) / (trisAtDepth+1);
  return textWidth;
}

float getHeightGap(int totalDepth, float tc, float bc) {
   float totalTreeHeight = WORLDHEIGHT - (tc + 2*textAscent() + bc);
   return totalTreeHeight / (totalDepth+1);
}

Point getTriPoint(MapTriangle mt, int totalDepth, float topCush, float botCush) {
  float heightGap = getHeightGap(totalDepth, topCush, botCush);
  float txtBoxWidth = getTextWidth(mt.depth);
  float offset = txtBoxWidth/2; // center it
  float xPos = TREEXSTART + offset + txtBoxWidth * mt.lvlRank;
  float yPos = WORLDHEIGHT - botCush - heightGap * mt.depth;
  return new Point(xPos, yPos);
}

Point getPolyPoint(int polyNum, float botCush) {
  float zoneStartx = TREEXSTART;
  float zoneStarty = WORLDHEIGHT - botCush;
  float zoneWidth = WORLDWIDTH - TREEXSTART;
  float zoneHeight = botCush - 2*textAscent();
  float slotSize = zoneWidth / (completedPolygons.size()+1); // don't forget outer face
  float offset = slotSize / 2;
  float xPos = zoneStartx + slotSize * (float) polyNum + offset;
  float yPos = zoneStarty + zoneHeight - textAscent();
  return new Point(xPos, yPos);
}

void drawTriTree() {
  float tc = 20; // distance from top of ceiling for the tree
  float bc = 75; // distance fromt he bottom of floor for the tree
  fill(255, 255, 0);
  rect(TREEXSTART, 0, WORLDWIDTH - TREEXSTART, WORLDHEIGHT);
  fill(0);
  strokeWeight(3);
  textAlign(CENTER);
  text("SUBDIVISION HIERARCHY", (WORLDWIDTH - TREEXSTART)/2 + TREEXSTART, 20);
  text("POLYGONS", (WORLDWIDTH - TREEXSTART)/2 + TREEXSTART, WORLDHEIGHT-textAscent());
  for (int i = 0; i < completedPolygons.size()+1; i++) {
    Point polyP = getPolyPoint(i, bc);
    if (result == i-1) fill(0, 255, 0);
    textAlign(CENTER);
    textSize(20);
    text(i-1, polyP.x, polyP.y);
    textSize(12);
    textAlign(0);
    fill(0);
  }
  textAlign(0);
  strokeWeight(1);
  for (MapTriangle mt : allTris) { // this will already know its rank and its depth
    Point p = getTriPoint(mt, treeDepth, tc, bc);
    fill(0);
    text(mt.id, p.x, p.y);
    if (mt.children.size() != 0) {
      for (MapTriangle child : mt.children) {
        Point childP = getTriPoint(child, treeDepth, tc, bc);
        if (QueryPoint != null || mt.depth == 1) stroke(0);
        else if (mt.depth == 2) stroke(255, 0, 0);
        else if (mt.depth == 3) stroke(0, 255, 0);
        else if (mt.depth == 4) stroke(128, 128, 0); // already using yellow
        else if (mt.depth == 5) stroke(0, 0, 255);
        else if (mt.depth == 6) stroke(255, 0, 255);
        else if (mt.depth == 7) stroke(0, 255, 255);
        else if (mt.depth == 8) stroke(255, 255, 255);
        if (mt == CURRTRI && child == NEXTTRI) {
          strokeWeight(10);
          stroke(255, 0, 0);
        }
        line(p.x, p.y, childP.x, childP.y - textAscent());
        strokeWeight(1);
      }
    }
    else { // make a line to each of their faceNums
      int polyNum = mt.faceNum;
      Point polyPoint = getPolyPoint(polyNum+1, bc);
      if (mt == CURRTRI && result < -1) {
        strokeWeight(10);
        stroke(0, 255, 0);
      }
      else stroke(0);
      textSize(20);
      line(p.x, p.y, polyPoint.x, polyPoint.y-textAscent());
      textSize(12);
      strokeWeight(1);
    }
  }
}

void writeToConsole(String message) {
  drawConsole();
  textSize(18);
  textAlign(CENTER);
  fill(0, 255, 0);
  text(message, TREEXSTART/2, CONSHEIGHT/2);
  fill(0);
  strokeWeight(1);
  textAlign(0);
  textSize(12);
}

void drawConsole() {
  stroke(255);
  strokeWeight(3);
  fill(0);
  rect(0, 0, TREEXSTART, CONSHEIGHT);
  fill(255);
  text("CONSOLE:", 5, 20);
  stroke(0);
  strokeWeight(1);
}

void removeIndPoint(HashMap<Point, HashMap<String, Object>> lookupTable) { // just removes triangles
  Iterator iter = indSet.iterator();
  Point p = (Point) iter.next();
  HashSet<MapTriangle> removedTris = removePoint(p, lookupTable);
  for (MapTriangle tri : removedTris) triangulation.remove(tri);
  ArrayList<Point> hole = getHole(p, lookupTable); // will be clockwise
  HashSet<MapTriangle> newTris = new HashSet<MapTriangle>();
  triangulatePolygon(hole, -2, newTris, true); // need to update adjacencies in lookuptable
  for (MapTriangle nt : newTris) {
    for (MapTriangle rem : removedTris) {
      if (trianglesOverlap(nt, rem)) {
        nt.addChild(rem);
      }
    }
    putTriInTree(nt);
    addTriToTable(nt, lookupTable);
    triangulation.add(nt);
  }
  indSet.remove(p);
}

void addToOuterPolygon(Point newP) { // TODO: just have 1 addToPolygon, do appropriate checks before
  if (!isInsideTriangle(newP, LEFT, TOP, RIGHT) || !noIntersects(newP, outerPolygon, completedPolygons)) return;
  outerPolygon.add(newP);
  idNum++;
}

void drawPoint(Point p) {
  if (p.equals(LEFT) || p.equals(RIGHT) || p.equals(TOP)) return;
  stroke(0);
  fill(0);
  ellipse(p.x, p.y, 10, 10);
}

void drawLineSegment(Point p1, Point p2) {
  stroke(0);
  line(p1.x, p1.y, p2.x, p2.y);
}

void drawLineSegment(Point p1, Point p2, int colNum) {
  stroke(colNum);
  line(p1.x, p1.y, p2.x, p2.y);
}

Point getSplitSpot(Point p) {
  Point res = snapToPoint(p, currentPolygon);
  if (res != null) return res;
  for (int i = 0; i < completedPolygons.size(); i++) {
    res = snapToPoint(p, completedPolygons.get(i));
    if (res != null) 
      return res;
  }
  return res;
}

void addPointToPolygon(Point p, ArrayList<Point> polygon) {
  polygon.add(p);
}

// big credits to Eliot Alter for this
void makeSplit(ArrayList<Point> split, int toSplitInd) {
  ArrayList<Point> toSplit = completedPolygons.get(toSplitInd);
  ArrayList<ArrayList<Point>> parts = new ArrayList<ArrayList<Point>>(); // parts = { {}, {} }
  parts.add(new ArrayList<Point>());
  parts.add(new ArrayList<Point>());
  int currPartInd = 0;
  ArrayList<Point> noEndPoints = new ArrayList<Point>();
  for (int ind = 1; ind < split.size()-1; ind++) 
    noEndPoints.add(split.get(ind)); // makes a copied list of just the middle vertices
  Point [] ends = {split.get(0), split.get(split.size()-1)};
  for (int i = 0; i < toSplit.size(); i++) { // iterate over the polygon
    if( toSplit.get(i).equals(ends[0]) ||
       toSplit.get(i).equals(ends[1]) ) { // if we find the intersection points
        parts.get(0).add(toSplit.get(i));
        parts.get(1).add(toSplit.get(i));
        ArrayList<Point> toAppend = noEndPoints;
        if (toSplit.get(i).equals(ends[1])) {
          toAppend = (ArrayList<Point>) noEndPoints.clone();
          Collections.reverse(toAppend);
        }
        for (int j = 0; j < toAppend.size(); j++) 
          parts.get(currPartInd).add(toAppend.get(j));
        currPartInd = (currPartInd == 0) ? 1 : 0;
     } else {
       parts.get(currPartInd).add(toSplit.get(i));
     }
  }
  completedPolygons.remove(toSplitInd);
  completedPolygons.add(parts.get(0));
  completedPolygons.add(parts.get(1));
}

void triangulate(ArrayList<ArrayList<Point>> polygons) {
  for (int i = 0; i < polygons.size(); i++) { //  TODO: clone the polygons so you can just remove points without messing up data
    int currFace = i;
    ArrayList<Point> polygon = polygons.get(i);
    triangulatePolygon(polygon, currFace, triangulation, isClockwise);
    for (MapTriangle mt : triangulation) {
      if (!allTris.contains(mt)) putTriInTree(mt);
    }
  }
  ArrayList<Point> convexHull = (ArrayList<Point>) outerPolygon.clone(); // after triangulatePockets...this turns into the convex hull
  triangulatePockets(convexHull, -1, triangulation, isClockwise);
  triangulateOutertriangle(convexHull, triangulation, isClockwise);
  for (MapTriangle mt : triangulation) {
    if (!allTris.contains(mt)) putTriInTree(mt);
  }
}

void triangulatePolygon(ArrayList<Point> polygon, int faceNum, HashSet<MapTriangle> tris, boolean clockwise) {
  // simple triangulation using ear clippings in O(n^3)... yuck
  if (polygon.size() < 3) {
    return; // can't be less than 3
  }
  int prevInd = polygon.size() - 1;
  int earInd = 0;
  int nextInd = 1;
  // find the first valid ear
  while (!canCutAtVertex(polygon, earInd, clockwise) ) { // because > 2 ears exist, this must terminate
    prevInd = (prevInd + 1) % polygon.size();
    earInd++; // won't ever wrap around
    nextInd = (nextInd + 1) % polygon.size();
  }
  Point [] triVerts = {polygon.get(prevInd), polygon.get(earInd), polygon.get( nextInd) };
  MapTriangle newT =  new MapTriangle(triVerts, faceNum, triNum);
  tris.add( newT); // make the triangulation
  triNum++;
  ArrayList<Point> withoutEar = new ArrayList<Point>(); // new polygon without the ear we just found
  for (int i = 0; i < polygon.size(); i++) { // maybe clone the polygon before this method
    Point p = polygon.get(i);
    if (!p.equals(polygon.get(earInd))) withoutEar.add(p);
  }
  triangulatePolygon(withoutEar, faceNum, tris, clockwise);
}

void triangulatePockets(ArrayList<Point> polygon, int faceNum, HashSet<MapTriangle> tris, boolean clockwise) { // find the "anti-ears" of the outerpolygon
  if (polygon.size() < 3) return;
  //boolean clockwise = orientationIsClockwise(polygon); // TODO: do this before hand and pass it in
  int prevInd = polygon.size() - 1;
  int currInd = 0;
  int nextInd = 1;
  // find the first anti-ear
  while (currInd < polygon.size() && !canCutAtVertex(polygon, currInd, !clockwise)) { // invert clockwise so it thinks outside is inside
    prevInd = (prevInd + 1) % polygon.size();
    currInd++; // won't ever wrap around
    nextInd = (nextInd + 1) % polygon.size();
  }
  if (currInd == polygon.size()) return; // there aren't anymore pockets... the polygon is convex
  // other wise we found a reflex vertex
  Point [] triVerts = {polygon.get(prevInd), polygon.get(currInd), polygon.get( nextInd) };
  MapTriangle newT = new MapTriangle(triVerts, faceNum, triNum);
  tris.add( newT ); // make the triangulation
  triNum++;
  polygon.remove(currInd);
  triangulatePockets(polygon,faceNum, tris, clockwise);
}

void triangulateOutertriangle(ArrayList<Point> convexHull, HashSet<MapTriangle> tris, boolean clockwise) {
  int minXind = getMinXindex(convexHull);
  int maxXind = getMaxXindex(convexHull);
  ArrayList<Point> outerPoly1 = new ArrayList<Point>();  
  ArrayList<Point> outerPoly2 = new ArrayList<Point>();
  outerPoly1.add(LEFT);
  outerPoly2.add(LEFT);
  int ind = minXind;
  while (ind != maxXind){
    outerPoly1.add(convexHull.get(ind));
    if (clockwise) ind = (ind+1) % convexHull.size();
    else ind = (ind == 0) ? convexHull.size()-1 : ind-1;
  }
  outerPoly1.add(convexHull.get(ind));
  outerPoly1.add(RIGHT);
  outerPoly1.add(TOP);
  ind = minXind;
  while (ind != maxXind) {
    outerPoly2.add(convexHull.get(ind));
    if (clockwise) ind = (ind == 0) ? convexHull.size()-1 : ind-1; 
    else ind = (ind+1) % convexHull.size();
  }
  outerPoly2.add(convexHull.get(ind));
  outerPoly2.add(RIGHT);
  triangulatePolygon(outerPoly1, -1, tris, false);
  triangulatePolygon(outerPoly2, -1, tris, true);
}

void makeTable(HashMap<Point, HashMap<String, Object>> table, HashSet<MapTriangle> tris) { // TODO: modify this
  for (MapTriangle tri : tris) {
    makeEdge(table, tri.verts[0], tri.verts[1]);
    makeEdge(table, tri.verts[1], tri.verts[2]);
    makeEdge(table, tri.verts[2], tri.verts[0]);
    addTriangle(table, tri.verts[0], tri); 
    addTriangle(table, tri.verts[1], tri);
    addTriangle(table, tri.verts[2], tri);
  }
}

void makeEdge(HashMap<Point, HashMap<String, Object>> table, Point p, Point q) {
  if (!table.containsKey(p)) {
    HashMap<String, Object> adjData = new HashMap<String, Object>();
    adjData.put("adjNodes", new HashSet<Point>());
    adjData.put("adjTris", new HashSet<MapTriangle>());
    table.put(p, adjData);
  }
  if (!table.containsKey(q)) {
    HashMap<String, Object> adjData = new HashMap<String, Object>();
    adjData.put("adjNodes", new HashSet<Point>());
    adjData.put("adjTris", new HashSet<MapTriangle>());
    table.put(q, adjData);
  }
  HashSet<Point> padjNodes = (HashSet<Point>) table.get(p).get("adjNodes");
  padjNodes.add(q);
  HashSet<Point> qadjNodes = (HashSet<Point>) table.get(q).get("adjNodes");
  qadjNodes.add(p);
}

void addTriangle(HashMap<Point, HashMap<String, Object>> table, Point p, MapTriangle tri) { // the triangle set must already exist
  HashSet<MapTriangle> adjTris = (HashSet<MapTriangle>) table.get(p).get("adjTris");
  adjTris.add(tri);
}

void getIndSet(HashMap<Point, HashMap<String, Object>> lookupTable, HashSet<Point> marked) {
  int degree = 8;
  Iterator it = lookupTable.entrySet().iterator();
  while (it.hasNext()) {
      Map.Entry pair = (Map.Entry)it.next();
      Point p = (Point) pair.getKey();
      if (!marked.contains(p)) {
        HashMap<String, Object> adj = (HashMap<String, Object>) pair.getValue();
        HashSet<Point> adjNodes = (HashSet<Point>) adj.get("adjNodes");
        if (adjNodes.size() <= degree) {
          indSet.add(p);
          fill(255, 0, 0);
          ellipse(p.x, p.y, 10, 10);
          for (Point q : adjNodes) marked.add(q);
        }
      }
  }
}

ArrayList<Point> getHole(Point p, HashMap<Point, HashMap<String, Object>> lookupTable) {
  HashSet<Point> padjNodes = (HashSet<Point>) lookupTable.get(p).get("adjNodes");
  ArrayList<Point> newPolygon = new ArrayList<Point>();
  for (Point q : padjNodes) newPolygon.add(q);
  radialSort(newPolygon, p); // sorts the points in counterclockwise order
  return newPolygon;
}

HashSet<MapTriangle> removePoint(Point p, HashMap<Point, HashMap<String, Object>> lookupTable) {
  HashSet<Point> adjNodes = (HashSet<Point>) lookupTable.get(p).get("adjNodes");
  HashSet<MapTriangle> adjTris = (HashSet<MapTriangle>) lookupTable.get(p).get("adjTris");
  for (Point neighbor : adjNodes) {
    HashSet<Point> neighborAdjNodes = (HashSet<Point>) lookupTable.get(neighbor).get("adjNodes");
    HashSet<MapTriangle> neighborAdjTris = (HashSet<MapTriangle>) lookupTable.get(neighbor).get("adjTris");
    neighborAdjNodes.remove(p);
    for (MapTriangle tri : adjTris) {
      if (neighborAdjTris.contains(tri)) neighborAdjTris.remove(tri);
    }
  }
  //lookupTable.remove(p); // delete the point from the table
  return adjTris;
}

void addTriToTable(MapTriangle mt, HashMap<Point, HashMap<String, Object>> lookupTable) {
  // we really just need one makeEdge
  makeEdge(lookupTable, mt.verts[0], mt.verts[1]);
  makeEdge(lookupTable, mt.verts[1], mt.verts[2]);
  makeEdge(lookupTable, mt.verts[2], mt.verts[0]);
  addTriangle(lookupTable, mt.verts[0], mt); 
  addTriangle(lookupTable, mt.verts[1], mt);
  addTriangle(lookupTable, mt.verts[2], mt);
}