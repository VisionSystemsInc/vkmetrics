#ifndef vkm_histogram_h
#define vkm_histogram_h

//:
// \file
// \brief A simple vkm_histogram class
// \author J.L.Mundy
// \date  June 3, 2017
//
// \verbatim
//  Modifications
//   <none yet>
// \endverbatim

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

template <class T>
class vkm_histogram{

public:

  vkm_histogram(){}
  vkm_histogram(size_t nbins, T min_val, T max_val);
  vkm_histogram(T min_val, T max_val, std::vector<T> const& counts){
    *this = vkm_histogram(counts.size(), min_val, max_val);
    counts_ = counts;
  }
  void upcount(T val, T weight = T(1));
  T bin_count(size_t bin) const{
    assert(bin<nbins_); return counts_[bin];}
  T count(T val) const;
  std::vector<T> count_array() const{return counts_;}
  T area() const;
  size_t nbins() const {return nbins_;}
  T min_val() const {return min_val_;}
  T max_val() const {return max_val_;}
  //: bin interval
  T delta() const {return delta_;}
  //: The average value for a bin
  T avg_bin_value(const size_t bin) const
  { assert(bin<nbins_); return min_val_ + bin*delta_ + delta_/2; }

  void print(std::ostream& os) const;

  //: Value for area fraction below value
  T value_with_area_below(const T area_fraction) const{
    if (area_fraction>T(1))
        return 0;
    T a = area();
    if (a == T(0))
      return 0;
    T sum = 0;
    for (unsigned int i=0; i<nbins_; ++i)
      {
        sum += counts_[i];
        if (sum>=area_fraction*a)
          return (i+1)*delta_+min_val_;
      }
    return 0;
  }
  //: Value for area fraction below value
  T value_with_area_above(const T area_fraction) const{
    if (area_fraction>T(1))
      return 0;
    T a = area();
    if (a == T(0))
      return 0;
    T sum = 0;
    for (unsigned int i=nbins_-1; i!=0; i--)
      {
        sum += counts_[i];
        if (sum>area_fraction*a)
          return (i+1)*delta_+min_val_;
      }
    return 0;
  }
  T most_probable_value() const{
    T max_v = T(0);
    size_t max_bin;
    for(size_t i = 0; i<nbins_; ++i)
      if(counts_[i]>max_v){
        max_v = counts_[i];
        max_bin = i;
      }
    return avg_bin_value(max_bin);
  }

private:

  size_t nbins_;
  T min_val_;
  T max_val_;
  T delta_;
  T range_;
  std::vector<T> counts_;

};

template <class T>
vkm_histogram<T>::vkm_histogram(size_t nbins, T min_val, T max_val):min_val_(min_val), max_val_(max_val){
  if (nbins>0)
  {
    nbins_ = nbins;
    range_ = max_val_-min_val_;
    delta_ = range_/nbins;
    counts_.resize(nbins, T(0));
  }
  else
  {
    nbins_ = 0;
    range_ = 0;
    delta_ = 0;
  }
}

template <class T>
void vkm_histogram<T>::upcount(T x, T mag)
{
  if (x<min_val_||x>max_val_)
    return;
  for (unsigned int i = 0; i<nbins_; ++i)
    if (T((i+1)*delta_) + min_val_ >= x)
    {
      counts_[i] += mag;
      break;
    }
}

template <class T>
T vkm_histogram<T>::count(const T val) const
{
  if (val<min_val_||val>max_val_)
    return T(0);
  for (size_t i = 0; i<nbins_; ++i)
    if (T((i+1)*delta_) + min_val_ >= val)
      return counts_[i];
  return T(0);
}

template <class T>
T vkm_histogram<T>::area() const{
  T a = T(0);
  for(size_t i = 0;  i<nbins_; ++i)
    a += counts_[i];
  return a;
}

template <class T>
void vkm_histogram<T>::print(std::ostream& os) const
{
  T a = this->area();
  if (a == T(0)) a = T(1);
  for (size_t i=0; i<nbins_; ++i)
    os << avg_bin_value(i) << ' ' << counts_[i]/a << std::endl;
}

#endif
