#ifndef vkm_containment_tree_h_
#define vkm_containment_tree_h_

#include <iostream>
#include <iosfwd>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <utility>
#include <algorithm>

// VXL ------------------------------------
#include <vcl_compiler.h>
#include <vgl/vgl_box_2d.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_polygon.h>
// ----------------------------------------


class cont_tree_node;// forward declaration

// the node of a region containment tree
class cont_tree_node{

public:

  cont_tree_node(size_t region_indx):region_indx_(region_indx){}
  void add_child(std::shared_ptr<cont_tree_node> const& parent, size_t region_indx){
    std::shared_ptr<cont_tree_node> temp(new cont_tree_node(region_indx));
    temp->parent_ = parent;
    child_nodes_.push_back(temp);
  }
  bool is_leaf(){
    return child_nodes_.size() == 0;
  }
  bool is_root(){
    return !parent_;
  }
  // a recursive function to print the tree
  static void print(std::shared_ptr<cont_tree_node> const & node){
    if(!node)
      return;
    if(node->child_nodes_.size()>0){
      std::cout << "self " << node->region_indx_ << ' ' << std::endl;
      std::cout << " children " ;
      for(std::vector<std::shared_ptr<cont_tree_node> >::const_iterator cit =  node->child_nodes_.begin();
          cit != node->child_nodes_.end(); ++cit)
        std:: cout << (*cit)->region_indx_ << ' ';
      std::cout << std::endl;
    }else return;
    for(std::vector<std::shared_ptr<cont_tree_node> >::const_iterator cit =  node->child_nodes_.begin();
        cit != node->child_nodes_.end(); ++cit)
      print(*cit); //recursive call on children
  }

  // members
  std::shared_ptr<cont_tree_node> parent_;                    // the parent tree node
  size_t region_indx_;                                        // the region index as read from input
  std::vector<std::shared_ptr<cont_tree_node> > child_nodes_; // the list of children

};

// a class to represent a multiply connected region
// the vgl_polygon can also store multiple sheets to represent holes
// but it is desirable to maintin the region index reference, thus the special class
// note that proper vertex order is *not* observed i.e. outer boundary counter clockwise,
// inner hole boundaries clockwise
class mc_region_2d{

public:

  mc_region_2d(): enclosing_region_(-1){}

  // convert mc_region_2d to a vgl_polygon to support pixel scanning
  vgl_polygon<double> poly() const{
    vgl_polygon<double> ret;
    ret.push_back(outer_cycle_);
    for(std::map<size_t, std::vector<vgl_point_2d<double> > >::const_iterator hit = holes_.begin();
        hit != holes_.end(); ++hit)
      ret.push_back(hit->second);
    return ret;
  }

  //:members
  // the outer boundary of the region
  std::vector<vgl_point_2d<double> > outer_cycle_;
  // the interior holes
  std::map<size_t, std::vector<vgl_point_2d<double> > > holes_;

  // the id of the enclosing region.
  // enclosing region not defined is indicated by -1
  size_t enclosing_region_;

};

std::ostream&  operator <<(std::ostream& os, mc_region_2d const& mcr_2d);

// the 3-d version of mc_region_2d
// all the vertices lie in the same 3-d plane
class mc_region_3d{

public:

  mc_region_3d():enclosing_region_(-1){}

  //: create a 2d multiply connected region from a projection of the 3d region onto the x-y plane
  mc_region_2d xy_region() const{
    std::vector<vgl_point_2d<double> > outer_cycle_2d;
    std::map<size_t, std::vector<vgl_point_2d<double> > > hole_map_2d;

    for(std::vector<vgl_point_3d<double> >::const_iterator v3dit = outer_cycle_.begin();
        v3dit != outer_cycle_.end(); ++v3dit)
      outer_cycle_2d.push_back(vgl_point_2d<double>(v3dit->x(), v3dit->y()));

    for(std::map<size_t, std::vector<vgl_point_3d<double> > >::const_iterator hit = holes_.begin();
        hit != holes_.end(); ++hit){
      size_t index = hit->first;
      const std::vector<vgl_point_3d<double> >& hverts = hit->second;

      std::vector<vgl_point_2d<double> > hole_verts_2d;
      for(std::vector<vgl_point_3d<double> >::const_iterator v3dit = hverts.begin();
          v3dit != hverts.end(); ++v3dit)
        hole_verts_2d.push_back(vgl_point_2d<double>(v3dit->x(), v3dit->y()));

      hole_map_2d[index] = hole_verts_2d;
    }
    mc_region_2d mcr_2d;
    mcr_2d.outer_cycle_ = outer_cycle_2d;
    mcr_2d.holes_ = hole_map_2d;
    mcr_2d.enclosing_region_ = enclosing_region_;
    return mcr_2d;
  }

  //: members
  std::vector<vgl_point_3d<double> > outer_cycle_;
  std::map<size_t, std::vector<vgl_point_3d<double> > > holes_;
  size_t enclosing_region_;

};

class vkm_containment_tree{

public:

  vkm_containment_tree(): verbose_(false){}

  // the polygons all have one sheet -- a simple boundary
  // but may have containment relations with other polys
  vkm_containment_tree(std::vector<vgl_polygon<double> > const& polys);

  //individual vertex sheets
  vkm_containment_tree(std::vector<std::vector<vgl_point_2d<double> > > sheets);

  // the polygons all have one sheet -- simple boundaries, relate to an index
  vkm_containment_tree(std::map<size_t, vgl_polygon<double> > const& poly_map);

  // the polygons all have one sheet -- simple boundaries
  // precomputed bounding box
  vkm_containment_tree(std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  > const& poly_map):
  input_regions_(poly_map){}
  void set_verbose(bool verbose){verbose_ = verbose;}

  //: form the containment maps by checking polygon inclusion
  void construct_containment_maps();

  //: set the roots of containment trees
  void set_root_nodes();
  void set_isolated_roots();
  void print_roots() {
	  for (std::vector<std::shared_ptr<cont_tree_node> >::iterator rit = roots_.begin();
		  rit != roots_.end(); ++rit) {
		  cont_tree_node::print(*rit);
	  }
  }

  //: a recursive method to form the region containment tree
  bool add_child(std::shared_ptr<cont_tree_node > & parent);

  //: build a region containment tree starting at a root node
  void build_containment_tree(std::shared_ptr<cont_tree_node >& root){
    while (add_child(root))
      continue;
  }
  void build_cont_trees();

  //: using the containment trees - construct mulitply connected (mc) regions
  //: all the points interior to a mc region are coplanar
  void mc_regs_recursive(std::shared_ptr<cont_tree_node >& node);

  //: using the containment tree form mc regions from a node and its immediate children
  void construct_multiply_connected_regions();

  //: main process method - carries out the
  //  formation of the containment tree and mc regions
  void process(){
    construct_containment_maps();
    set_root_nodes();
    set_isolated_roots();
    build_cont_trees();
    construct_multiply_connected_regions();
  }
  // the multiply connected regions formed from the containment tree
  std::map < size_t,  mc_region_2d>& mc_regions(){ return  multiply_connected_regs_;}

private:

  bool verbose_;// print debug output

  std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  > input_regions_;

  //: regions that improperly intersect
  std::vector<size_t> intersecting_regions_;

  //: the "contains" relation between regions
  std::map<size_t, std::vector<size_t> > contains_;

  //: the "contained by" relation between regions
  std::map<size_t, std::vector<size_t> > contained_by_;

  //: regions neither contained or contained by
  std::vector<size_t> isolated_regions_;

  //: the roots of the region containment trees - isolated regions will just have a root node
  std::vector<std::shared_ptr<cont_tree_node> > roots_;

 //: the mulitply connected regions in image coordinates
 std::map < size_t,  mc_region_2d> multiply_connected_regs_;

};

#endif
