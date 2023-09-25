
#include <map>
#include <iostream>

typedef std::map<int, int> map_type;

extern "C"{
  void c_map_init(map_type** p){
    *p = new map_type();
  }

  void c_map_insert(map_type** p, int* k, int* v){
    (*p)->insert(std::pair<int,int>(*k, *v));
  }

  void c_map_final(map_type** sd){
    if(*sd != nullptr){
      std::map<int, int>().swap(**sd);
      *sd = nullptr;
    }
  }

  void c_map_size(map_type** p, int* v){
    *v = (*p)->size();
  }

  void c_map_find(map_type** p, int* key, int* v){
    map_type::iterator it;
    it = (*p)->find(*key);
    
    if(it != (*p)->end()){
      *v = it->second;
    }else{
      *v = -1;
    }
  }
  
}
