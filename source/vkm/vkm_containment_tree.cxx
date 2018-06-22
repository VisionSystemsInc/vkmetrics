#include <math.h>
#include <cassert>

// VXL ------------------------------------
#include <vgl/vgl_area.h>
#include <vgl/vgl_distance.h>
#include <vgl/vgl_intersection.h>
#include <vgl/vgl_clip.h>
// ----------------------------------------

#include "vkm_containment_tree.h"


std::ostream& operator <<(std::ostream& os, mc_region_2d const& mcr_2d)
{
  vgl_polygon<double> ply = mcr_2d.poly();//looses original hole id FIXME
  os << ply;
  return os;
}


vkm_containment_tree::vkm_containment_tree(
    std::vector<vgl_polygon<double> > const& polys)
{
  std::map<size_t, vgl_polygon<double> > poly_map;
  for(size_t i = 0; i<polys.size(); ++i)
    poly_map[i] = polys[i];
  *this = vkm_containment_tree(poly_map);
}


vkm_containment_tree::vkm_containment_tree(
    std::map<size_t, vgl_polygon<double> > const& poly_map)
{
  for (std::map<size_t, vgl_polygon<double> >::const_iterator pit = poly_map.begin();
       pit != poly_map.end(); ++pit) {
    size_t indx = pit->first;
    size_t n_sheets = pit->second.num_sheets();
    assert(n_sheets == 1);
    const std::vector<vgl_point_2d<double> >& sheet = pit->second[0];
    size_t n = sheet.size();
    vgl_box_2d<double> bb;
    for (size_t i = 0; i < n; ++i)
      bb.add(sheet[i]);
    std::pair<vgl_polygon<double>, vgl_box_2d<double> > pr(pit->second, bb);
    input_regions_[indx] = pr;
  }
}


//for use with 2D Hull refinement (improperly in modl_hull_facade_2d for now)
vkm_containment_tree::vkm_containment_tree(
    std::vector<std::vector<vgl_point_2d<double> > > sheets)
{
  for(size_t i = 0; i<sheets.size();i++)
  {
    std::vector<vgl_point_2d<double> > currSheet = sheets.at(i);
    vgl_box_2d<double> bb;
    for(size_t j = 0; j < currSheet.size();j++)
    {
      bb.add(currSheet.at(j));
    }
    vgl_polygon<double> currPoly(currSheet);
    if(vgl_area(bb) == 0)
    {
      std::cout << "zero area region[" << i << "] - fatal " << std::endl;
    }
    else
    {
      std::pair<vgl_polygon<double>, vgl_box_2d<double> > pr(currPoly, bb);
      //should be renamed because this doesnt apply
      input_regions_[i] = pr;
    }
  }
}


static bool valid_overlap(
    double overlap_tol, vgl_polygon<double> const& poly_i,
    vgl_polygon<double> const& poly_j)
{
  vgl_polygon<double> int_poly = vgl_clip(poly_i, poly_j);
  // find longest edge
  size_t n_sheets = int_poly.num_sheets();
  if(n_sheets == 0)
    return true;
  // find most distant vertices in the intersection
  // typically the length of a long thin overlap region
  double max_dist = 0.0;
  for(size_t j = 0; j<n_sheets; ++j){
    for(size_t i = 0; i < int_poly[j].size(); ++i){
      const vgl_point_2d<double>& pi = int_poly[j][i];
      for(size_t k = 0; k<n_sheets; ++k){
        for(size_t m = 0; m < int_poly[k].size(); ++m){
          if(j==k && m<=i)
            continue;
          const vgl_point_2d<double>& pm = int_poly[k][m];
          vgl_vector_2d<double> v = pi-pm;
          double len = v.length();
          if(len > max_dist)
            max_dist = len;
        }
      }
    }
  }
  if(max_dist == 0.0)
    return true;
  // long thin intersection regions are allowed
  // just due to slight boundary overlap
  double area = vgl_area(int_poly);
  if(area/max_dist > overlap_tol)
    return false;
  return true;
}


void vkm_containment_tree::construct_containment_maps()
{
  for (std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  >::iterator riti = input_regions_.begin();
       riti != input_regions_.end(); ++riti) {
    size_t indx_i = riti->first;
    bool isolated = true;
    const std::pair<vgl_polygon<double>, vgl_box_2d<double> >& pr_i = riti->second;
    const vgl_polygon<double> & poly_i = pr_i.first;
    const vgl_box_2d<double>& bb_i = pr_i.second;
    for (std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  >::iterator ritj = input_regions_.begin();
         ritj != input_regions_.end(); ++ritj) {
      size_t indx_j = ritj->first;
      if (indx_j <= indx_i)// assume map sort is accending on indx
        continue;
      const std::pair<vgl_polygon<double>, vgl_box_2d<double> >& pr_j = ritj->second;
      const vgl_polygon<double>& poly_j = pr_j.first;
      const vgl_box_2d<double>& bb_j = pr_j.second;
      vgl_box_2d<double> bint = vgl_intersection(bb_i, bb_j);
      if (bint.is_empty())
        continue;

      // bounding boxes intersect so check if poly_j is inside poly_i
      // or if poly_i is inside poly_j
      const std::vector<vgl_point_2d<double> >& verts_i = poly_i[0];
      const std::vector<vgl_point_2d<double> >& verts_j = poly_j[0];

      bool all_contained = true;
      for (std::vector<vgl_point_2d<double> >::const_iterator vit = verts_j.begin();
           vit != verts_j.end() && all_contained; ++vit)
        if (!poly_i.contains(*vit))
          all_contained = false;

      bool j_contained_in_i = all_contained;

      all_contained = true;
      for (std::vector<vgl_point_2d<double> >::const_iterator vit = verts_i.begin();
           vit != verts_i.end() && all_contained; ++vit)
        if (!poly_j.contains(*vit))
          all_contained = false;
      bool i_contained_in_j = all_contained;

        // normal intersection
        if (j_contained_in_i){
          contains_[indx_i].push_back(indx_j);
          contained_by_[indx_j].push_back(indx_i);
          isolated = false;
        }
        if (i_contained_in_j){
          contains_[indx_j].push_back(indx_i);
          contained_by_[indx_i].push_back(indx_j);
          isolated = false;
        }
        if (indx_i == 16 && indx_j == 17)
          std::cout << ' ';
        // if the boundaries potentially intersect
        if (!i_contained_in_j && !j_contained_in_i) {
          if(valid_overlap(2.0,poly_i, poly_j)){
            isolated = true;
            continue;
          }else{
            //remove the smaller region
            double area_i = vgl_area(poly_i), area_j = vgl_area(poly_j);
            if(area_i > area_j){
              std::vector<size_t>::iterator iit;
              iit = std::find(intersecting_regions_.begin(), intersecting_regions_.end(), indx_j);
              if(iit == intersecting_regions_.end())
                  intersecting_regions_.push_back(indx_j);
            }else{
              std::vector<size_t>::iterator iit;
              iit = std::find(intersecting_regions_.begin(), intersecting_regions_.end(), indx_i);
              if (iit == intersecting_regions_.end())
                intersecting_regions_.push_back(indx_i);
              isolated = false;
            }
            std::cout << "two regions," << indx_i << " and " << indx_j << " intersect but not properly contained -- warning" << std::endl;
            continue;
            }
        }
    }//end indx_j
    if (isolated)
      isolated_regions_.push_back(indx_i);
  }// end indx_i

  // remove intersecting regions from maps
  for(std::vector<size_t>::iterator iit = intersecting_regions_.begin();
      iit != intersecting_regions_.end(); ++iit){
      contains_.erase(*iit);
      contained_by_.erase(*iit);
  }

  for(std::map<size_t, std::vector<size_t> >::iterator cit = contains_.begin();
      cit != contains_.end(); ++cit){
    size_t index = cit->first;
    std::vector<size_t> temp;
    const std::vector<size_t>& children= cit->second;
    for(std::vector<size_t>::const_iterator chit = children.begin();
        chit != children.end(); ++chit){
      std::vector<size_t>::iterator rit;
      rit = std::find(intersecting_regions_.begin(), intersecting_regions_.end(), *chit);
      if(rit == intersecting_regions_.end())
        temp.push_back(*chit);
    }
    cit->second = temp;
  }

  // remove any regions contained by the region being removed
  std::cout << "==containment ==" << std::endl;
  for (std::map<size_t, std::vector<size_t> >::iterator cit = contains_.begin();
       cit != contains_.end(); ++cit) {
    std::cout << cit->first << "  contains => ";
    const std::vector<size_t>& inside = cit->second;
    if(inside.size() == 0)
      std::cout << "null";
    else
      for (std::vector<size_t>::const_iterator iit = inside.begin();
           iit != inside.end(); ++iit)
        std::cout << (*iit) << ' ';
    std::cout << std::endl;
  }

  std::cout << "==being contained ==" << std::endl;
  for (std::map<size_t, std::vector<size_t> >::iterator cit = contained_by_.begin();
       cit != contained_by_.end(); ++cit) {
    std::cout << cit->first << " is contained by => ";
    const std::vector<size_t>& inside = cit->second;
    if (inside.size() == 0)
      std::cout << "null";
    else
      for (std::vector<size_t>::const_iterator iit = inside.begin();
           iit != inside.end(); ++iit)
        std::cout << (*iit) << ' ';
    std::cout << std::endl;
  }

}


void vkm_containment_tree::set_root_nodes()
{
  for (std::map<size_t, std::vector<size_t> >::iterator cit = contains_.begin();
       cit != contains_.end(); ++cit) {
    //is region contained by any other region?
    const std::vector<size_t>& cont_by = contained_by_[cit->first];
    if(cont_by.size() == 0){
      // is a root of a containment tree
      std::shared_ptr<cont_tree_node> root(new cont_tree_node(cit->first));
      roots_.push_back(root);
    }
  }
}


void vkm_containment_tree::set_isolated_roots()
{
	for (std::vector<size_t>::iterator iit = isolated_regions_.begin();
		iit != isolated_regions_.end(); ++iit) {
		if (contained_by_[*iit].size() != 0)
			continue;
		std::shared_ptr<cont_tree_node> iso_node(new cont_tree_node(*iit));
		roots_.push_back(iso_node);
	}
}


bool vkm_containment_tree::add_child(
    std::shared_ptr<cont_tree_node > & parent)
{
  //contains_ <<>> the set of regions that contain at least one region (not self)
  const std::vector<size_t>& cntnd = contains_[parent->region_indx_];
  if(cntnd.size() == 0)
    return false;
  for(std::vector<size_t>::const_iterator cit = cntnd.begin();
      cit != cntnd.end(); ++cit){
    bool found = false;

    // see if *cit is already a child of this parent
    for(std::vector<std::shared_ptr<cont_tree_node> >::iterator ecit = parent->child_nodes_.begin();
        ecit != parent->child_nodes_.end()&&!found; ++ecit) {
      if( (*ecit)->region_indx_ == *cit)
        found = true;
    }
	  if (found)
      continue;

    // parent contains *cit but it is not already a child
    // next check if *cit is contained by parent
    const std::vector<size_t>& cntnd_by = contained_by_[*cit];
    if (cntnd_by.size() == 0)
      continue; //no parent
    std::vector<size_t>::const_iterator bit;
    bit = std::find(cntnd_by.begin(), cntnd_by.end(), parent->region_indx_);
    if (bit == cntnd_by.end())
      continue; // parent doesn't contain *cit

    //has correct parent
    // now see if *cit is already contained by an existing child at this level
    bool contains = false;
    for (std::vector<std::shared_ptr<cont_tree_node> >::iterator ecit = parent->child_nodes_.begin();
         ecit != parent->child_nodes_.end() && !contains; ++ecit) {
      size_t ri = (*ecit)->region_indx_; // child index
      const std::vector<size_t>& cntns = contains_[ri]; // regions contained by child
      bit = std::find(cntns.begin(), cntns.end(), *cit);
      contains = bit != cntns.end();
    }
    if (contains) // another child already contains (*cit)
      continue;

    // all ok we can add child
    std::shared_ptr<cont_tree_node> temp( new cont_tree_node(*cit));
    temp->parent_ = parent;// set the parent link

    //recursive call - tree is constructed depth first
    while (add_child(temp))
      continue;

    // after recursive return, add child to parent node
    parent->child_nodes_.push_back(temp);
    return true;
  }

  return false;
}


void vkm_containment_tree::build_cont_trees()
{
	size_t n = roots_.size();
	for (size_t i = 0; i < n; ++i) {
    build_containment_tree(roots_[i]);
    std::cout << std::endl << "Tree: " << roots_[i]->region_indx_ << std::endl;
    cont_tree_node::print(roots_[i]);
  }
}


void vkm_containment_tree::mc_regs_recursive(
    std::shared_ptr<cont_tree_node >& node)
{
  size_t indx = node->region_indx_;

  size_t parent_index = -1;
  if(!node->is_root())
    parent_index = node->parent_->region_indx_;

  // case I no children
  if(node->is_leaf()){
    mc_region_2d leaf_mcr;
    leaf_mcr.enclosing_region_ = parent_index;
    leaf_mcr.outer_cycle_ = input_regions_[indx].first[0];
    multiply_connected_regs_[indx] = leaf_mcr;
    return;
  }
  mc_region_2d mcr;
  mcr.outer_cycle_ = input_regions_[indx].first[0];
  mcr.enclosing_region_ = parent_index;

  std::vector<std::shared_ptr<cont_tree_node> > cnodes = node->child_nodes_;
  size_t n_childrn = cnodes.size();
  for(size_t i = 0; i<n_childrn; ++i){
    const std::vector<vgl_point_2d<double> >& hole = input_regions_[cnodes[i]->region_indx_].first[0];
    mcr.holes_[cnodes[i]->region_indx_] = hole;
  }
  multiply_connected_regs_[indx] = mcr;
  for(size_t i = 0; i<n_childrn; ++i)
    mc_regs_recursive(cnodes[i]);
  return;
}


void vkm_containment_tree::construct_multiply_connected_regions()
{
  size_t n = roots_.size();
  for (size_t i = 0; i < n; ++i)
    mc_regs_recursive(roots_[i]);
}

