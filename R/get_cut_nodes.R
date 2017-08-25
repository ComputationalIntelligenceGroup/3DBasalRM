#' Check terminal nodes to select those that have been cut and return the id of that nodes
#'
#' The function computes the 3D convex hull of the terminal nodes. Then search the normal faces that points to the same direction.
#' If there are several coplanar faces there is a cut and all that terminal nodes are delimiting a cut.
#'
#' @param neuron is the reconstruction in plain text with JSON format
#' @param n_points_plane is an integer denoting the minimum number of points in a plane to be considered a cutting plane
#'
#' @return an array with the id of the terminal nodes that belong to the tips of the cut dendrites
#'
#' @export
get_cut_nodes<-function(neuron, n_points_plane=5)
{
  #Get all the terminal nodes of a neuron
  terminal_nodes <- c_get_terminal_nodes(neuron$plain)

  #Get the coordinates and the id
  terminal_nodes_coords <- terminal_nodes[, 1:3]
  terminal_nodes_id <- terminal_nodes[, 4]

  #Compute convex hull of the terminal nodes
  convex_hull <- list()
  convex_hull$vb <- t(cbind(terminal_nodes_coords, 1))
  convex_hull$it <- t(convhulln(terminal_nodes_coords))
  class(convex_hull) <- "mesh3d"

  #Orientante all the face normals to point out the mesh
  convex_hull <- vcgClean(convex_hull, 7)
  face_normals <- t(facenormals(convex_hull)$normals)

  #Cluster face normals to find cut planes
  cl <- dbscan(face_normals, 0.3, minPts=1)

  #If there are at least 5 faces pointing in the same direction it is a cut plane
  cut_face <- which(table(cl$cluster) > n_points_plane)
  idx_cut_faces <- which(cl$cluster %in% cut_face)

  #Get the closest position to the convex hull for all the terminal nodes
  node_projection <- vcgClostKD(terminal_nodes_coords, convex_hull)

  #Those nodes that are not in the convex hull are assigned to the closest face
  nodes_inside_convex_hull <- which(node_projection$quality!=0)
  close_face_nodes_inside <- node_projection$faceptr[nodes_inside_convex_hull]

  #Those nodes inside the convex hull whose closest face is a cut face are regarded as cut
  cut_nodes_inside <- nodes_inside_convex_hull[which(close_face_nodes_inside %in% idx_cut_faces)]
  cut_nodes_convex_hull <- which(unlist(lapply(vcgVFadj(convex_hull), function(x) {any(x %in% idx_cut_faces)})))

  #Join cut nodes in the convex hull and inside the convex hull
  cut_nodes <- c(cut_nodes_convex_hull, cut_nodes_inside)

  return(terminal_nodes_id[cut_nodes])
}
