#ifndef Pivot_hpp
#define Pivot_hpp

// Cole Foster
// 2023-03-20
// Struct for holding pivot information for recursive definition

struct Pivot {
   public:
    Pivot(){};
    Pivot(unsigned int index, float const radius) : _index(index), _radius(radius){};
    ~Pivot(){};

    inline void addChild(Pivot& childPivot, float const distance) {
        _pivotDomain.push_back(childPivot);
        if (distance > _maxChildDistance) {
            _maxChildDistance = distance;
        }
    }
    inline void addChild(unsigned int childIndex, float const distance) {
        _childCount++;
        _pivotDomain.push_back(Pivot(childIndex, 0.0f));  // initialize point as a pivot with 0 radius
        if (distance > _maxChildDistance) {
            _maxChildDistance = distance;
        }
    }

    unsigned int _index = 0;
    float _radius = 0.0f;
    int _childCount = 0;
    float _maxChildDistance = 0.0f;
    std::vector<Pivot> _pivotDomain{};
};

#endif  // Pivot_hpp