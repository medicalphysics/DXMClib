

#include <cstdint>
#include <iostream>
#include <vector>
struct Node {

    // bits 0..1 : splitting dimension
    // bits 2..30 : offset bits
    // bit 31 (sign) : flag whether node is a leaf

    struct {
        std::uint32_t dim : 2; // dimensjon for branch
        std::uint32_t offset : 29; // offset to first child (branch) or to first element (leaf)
        std::uint32_t flag : 1; // Node is leaf (1) or branch (0)
    } dim_offset_flag;

    union {
        float split = 0; // Union split is for branches
        std::uint32_t nelements; // Union nelements is for number of items
    } split_nelements;

    Node()
    {
        dim_offset_flag.dim = 0;
        dim_offset_flag.offset = 0;
        dim_offset_flag.flag = 0;
    }

    std::uint32_t dim()
    {
        return dim_offset_flag.dim;
    }
    void setDim(std::uint32_t dim)
    {
        dim_offset_flag.dim = dim;
    }

    void setLeaf()
    {
        dim_offset_flag.flag = std::uint32_t { 1 };
    }
    std::uint32_t isLeaf()
    {
        return dim_offset_flag.flag;
    }

    void setOffset(std::uint32_t offset)
    {
        dim_offset_flag.offset = offset;
    }
    std::uint32_t offset()
    {
        return dim_offset_flag.offset;
    }
};

template <typename U>
class FlatKDTreeTest {
public:
    FlatKDTreeTest(const std::vector<U>& elements)
    {
        m_items = elements;
    }

    void build(int max_depth = 8)
    {
        std::vector<std::uint32_t> indices(m_items.size());
        std::iota(indices.begin(), indices.end());
        build(indices, max_depth);
    }

protected:
    void build(std::vector<std::uint32_t>& indices, int max_depth = 8)
    {
        std::vector<std::uint32_t> left;
        std::vector<std::uint32_t> right;
        left.reserve(indices.size());
        right.reserve(indices.size());

        m_nodes.reserve(indices.size());
        m_nodes.push_back({});
        m_nodes[0].setOffset(0);
        std::size_t currentNodeIdx = 0;

        while (max_depth > 0) {
            // handle current node
            auto split_dim = split_dimension(indices);
            auto split_val = split_value(indices, split_dim);
            for (auto idx : indices) {
                const auto& item = m_items[idx];
                // assigning side
                if (item[split_dim] < split_val)
                    left.push_back(idx);
                else
                    right.push_back(idx);
            }



            // update current node
            auto& cnode = m_nodes[currentNodeIdx];

            cnode.setDim(split_dim);
            cnode.set

                -- max_depth;
        }
    }

    std::uint32_t split_dimension(const std::vector<std::uint32_t>& indices)
    {
        std::array<double, 6> extent = {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };
        for (auto idx : indices) {
            const auto& item = m_items[idx];
            for (int i = 0; i < 3; ++i) {
                extent[i] = std::min(extent[i], item[i]);
                extent[i + 3] = std::max(extent[i + 3], item[i]);
            }
        }
        std::array<double, 3> s;
        for (int i = 0; i < 3; ++i) {
            s[i] = extent[i + 3] - extent[i];
        }
        if (s[0] > s[1]) {
            if (s[0] > s[2])
                return 0;
            else
                return 1;
        } else {
            if (s[1] > s[2])
                return 1;
            else
                return 2;
        }
    }

    float split_value(const std::vector<std::uint32_t>& indices, std::uint32_t dimension)
    {
        std::vector<float> s(indices.size());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            const auto idx = indices[i];
            const auto& item = m_items[idx];
            s[i] = static_cast<float>(item[dimension]);
        }
        std::sort(s.begin(), s.end());
        return s[s.size() / 2];
    }

    int figureOfMerit(const std::vector<std::uint32_t>& indices, std::uint32_t dimension, float value)
    {
        int fom = 0;
        int shared = 0;
        for (auto idx:indices) {
            const auto& item = m_items[idx];
            int side = 1;
            if (item[dimension] > value )
            const auto side = planeSide(item, planesep, dim);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }

private:
    std::vector<std::uint32_t> m_indices;
    std::vector<U> m_items;
    std::vector<Node> m_nodes;
};

int main()
{

    std::cout << "sizeof Node: " << sizeof(Node) << std::endl;

    std::vector<Node> nodes(1);

    auto& root = nodes.at(0);

    Node n;
    n.setDim(2);

    /*
    //n.setLeaf();
    std::cout << "Is leaf: " << n.isLeaf() << std::endl;

    std::cout << std::hex;

    std::cout << n.dim_offset_flag << std::endl;
    std::cout << std::uint32_t { 1 } << 31 << std::endl;
*/
    return 0;
}
