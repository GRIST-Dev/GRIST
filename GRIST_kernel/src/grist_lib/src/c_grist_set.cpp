
#include <set>
#include <iostream>
#include <algorithm>
#include <vector>

class set_vec{
    std::set<int> p;
    std::vector<int> v;
    
public:
    void insert(int i){
        if(p.find(i) == p.end()){
            p.insert(i);
            v.push_back(i);
        }
    }
    
    void final(){
        p.clear();
        std::vector<int>().swap(v);
    }

    int size(){
        return v.size();
    }
    
    bool find(int i){
        if(p.find(i) == p.end()){
            return false;
        }else{
            return true;
        }
    }
    
    void dump(int* dv){
        std::copy(v.begin(), v.end(), dv);
    }
    
    void sort(int* dv){
        std::sort(v.begin(), v.end());
        std::copy(v.begin(), v.end(), dv);
    }
};

extern "C"{
  void c_set_init(set_vec** p){
    *p = new set_vec();
  }

  void c_set_insert(set_vec** p, int* v){
    (*p)->insert(*v);
  }

  void c_set_final(set_vec** sd){
    if(*sd != nullptr){
      (*sd)->final();
      *sd = nullptr;
    }
  }

  void c_set_size(set_vec** p, int* v){
    *v = (*p)->size();
  }

  void c_set_find(set_vec** p, int* v, int* r){
    if((*p)->find(*v)){
      *r = 1;
    }else{
      *r = 0;
    }
  }
  
  void c_set_dump(set_vec** p, int* v){
      (*p)->dump(v);
  }

  void c_set_sort(set_vec** p, int* v){
      (*p)->sort(v);
  }
}
