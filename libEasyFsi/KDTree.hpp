#pragma once
/*-------------------------------------------------------------------------
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from
the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not
   be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution.
-------------------------------------------------------------------------*/

//!-------------------------------------------------------------
//! @file       KDTree.hpp
//!             The definition KDTree class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <algorithm>
#include <numeric>
#include <array>
#include <vector>
#include <omp.h>

#include "Assert.hpp"

namespace EasyLib {

    //! @brief A fast k-dimension tree for nearest neighbor searching based on ANN kd-tree.
    //! @tparam CoordType  The coordinate type, default is double.
    //! @tparam NDIM       The dimension of data, default is 3.
    //! @tparam SizeType   Data type used for indexing, default is size_t.
    template<typename CoordType = double, size_t NDIM = 3, typename SizeType = size_t>
    class KDTree
    {
    public:
        static constexpr size_t ND = NDIM;
        static constexpr size_t LO = 0;
        static constexpr size_t HI = 1;

        using coord_type = CoordType;
        using size_type = SizeType;
        using point_type = std::array<coord_type, ND>;

        static constexpr coord_type ERR = 0.001;
        static constexpr coord_type INF = std::numeric_limits<coord_type>::max();

        //! @brief default constructor
        KDTree() = default;

        //! @brief copy constructor
        KDTree(const KDTree& tree)
        {
            copy_from_(tree);
        }

        //! @brief move constructor
        KDTree(KDTree&& tree)noexcept
        {
            move_(std::move(tree));
        }

        //! @brief Constructor
        //! @param coords Point coordinates array£º[x0,y0,..., x1,y1,..., ..., xn,yn,...]
        //! @param npts   Point number
        //! @param persistent Is the coordinates array persistent storage? A copy will be maked if is false.
        //! @note Coincide points are not allowed!
        KDTree(const coord_type* coords, size_type npts, bool persistent = true)
            :KDTree()
        {
            create(coords, npts, persistent);
        }

        KDTree& operator = (const KDTree& tree)
        {
            copy_from_(tree);
            return *this;
        }

        KDTree& operator = (KDTree&& tree)noexcept
        {
            if (&tree == this)return *this;
            move_(std::move(tree));
            return *this;
        }

        //! @brief Destructor
        ~KDTree() = default;

        //! @brief Clear all data.
        void clear()
        {
            points_ = nullptr;
            root_ = nullptr;
            bnd_box_lo_.fill(0);
            bnd_box_hi_.fill(0);
            pidx_.clear();
            nodes_.clear();
            points_stg_.clear();
        }

        //! @brief Create the KDTree
        //! @param coords     Point coordinates array: [x0,y0,..., x1,y1,..., ..., xn,yn,...]
        //! @param npts       Point number
        //! @param persistent Is the coordinates array persistent storage? A copy will be maked if is false.
        //! @note Coincide points are not allowed!
        void create(const coord_type* coords, size_type npts, bool persistent = true)
        {
            // initialize
            init_(coords, npts, persistent);

            if (npts == 0)return;

            // construct bounding box
            bounding_box_();

            // build KDtree
            point_type bnd_box_lo = bnd_box_lo_;
            point_type bnd_box_hi = bnd_box_hi_;
            root_ = build_(pidx_.data(), npts, bnd_box_lo, bnd_box_hi);
        }

        //! @brief Finding nearest neighbors for one point
        //! @param [in]  q       Coordinates of query point.
        //! @param [in]  n_neigh Number of nearest neighbors
        //! @param [out] idx     Nearest neighbor list
        //! @param [out] dist_sq Nearest distance square
        //! @param [in]  eps     Tolerance for nearest neighbor searching. Default is ZERO, which means precision nearest.
        //! @return Number of actual nearest neighbor.
        size_type search(const coord_type* q, size_type n_neigh, size_type* idx, coord_type* dist_sq, coord_type eps = 0)const
        {
            ASSERT(n_neigh <= static_cast<size_type>(pidx_.size()));
            n_neigh = std::min(n_neigh, static_cast<size_type>(pidx_.size()));
            if (n_neigh <= 0)return 0;

            // query data structure
            KDQuery Q{
                q,                         // q
                n_neigh,                   // n_neigh
                0,                         // n_found
                dist_sq,                   // key
                idx,                       // info
                (1.0 + eps) * (1.0 + eps), // max_err
                0                          // n_visited
            };

            // compute distance from point to box
            coord_type box_dist = 0;
            for (size_t d = 0; d < ND; ++d) {
                if (q[d] < bnd_box_lo_[d]) {
                    auto delt = bnd_box_lo_[d] - q[d];
                    box_dist += delt * delt;
                }
                else if (q[d] > bnd_box_hi_[d]) {
                    auto delt = bnd_box_hi_[d] - q[d];
                    box_dist += delt * delt;
                }
            }

            search_(root_, Q, box_dist);
            return Q.n_found;
        }

        //! @brief Finding nearest neighbors for given points.
        //! @param [in]  n_pnts  Number of query points
        //! @param [in]  pnts    Coordinate array of query points: [x0,y0,..., x1,y1,..., ..., xn,yn,...]
        //! @param [in]  n_neigh Number of nearest neighbors of each query point
        //! @param [out] idx     Nearest neighbor list: [n0_0,n0_1,..., n1_0,n1_1,..., ..., nn_0,nn_1,...]
        //! @param [out] dist_sq Nearest distance square: [d0_0,d0_1,..., d1_0,d1_1,..., ..., dn_0,dn_1,...]
        //! @param [in]  eps     Tolerance for nearest neighbor searching. Default is ZERO, which means precision nearest.
        void search(size_type n_pnts, const coord_type* pnts, size_type n_neigh, size_type* idxs, coord_type* dist_sq, coord_type eps = 0)const
        {
            ASSERT(n_neigh <= static_cast<size_type>(pidx_.size()));
            n_neigh = std::min(n_neigh, static_cast<size_type>(pidx_.size()));
            if (n_pnts <= 0 || n_neigh < 0)return;

            KDQuery Q{
                    pnts,                        // q
                    n_neigh,                     // n_neigh
                    0,                           // n_found
                    dist_sq,                     // key
                    idxs,                        // info
                    (1.0 + eps) * (1.0 + eps),   // max_err
                    0                            // n_visited
            };

            for (size_type i = 0; i < n_pnts; ++i) {

                // compute distance from point to box
                coord_type box_dist = 0;
                for (size_t d = 0; d < ND; ++d) {
                    if (Q.q[d] < bnd_box_lo_[d]) {
                        auto delt = bnd_box_lo_[d] - Q.q[d];
                        box_dist += delt * delt;
                    }
                    else if (Q.q[d] > bnd_box_hi_[d]) {
                        auto delt = bnd_box_hi_[d] - Q.q[d];
                        box_dist += delt * delt;
                    }
                }
                Q.n_found = 0;

                // search
                search_(root_, Q, box_dist);

                // fill other
                for (size_type j = Q.n_found; j < Q.n_neigh; ++j) {
                    Q.key[j] = INF;
                    Q.info[j] = -1;
                }

                Q.q += ND;
                Q.n_found = 0;
                Q.key += n_neigh;
                Q.info += n_neigh;
                Q.n_visited = 0;
            }
        }

        //! @brief Using OpenMP to finding nearest neighbors.
        //! @param [in]  n_pnts  Number of query points
        //! @param [in]  pnts    Coordinate array of query points: [x0,y0,..., x1,y1,..., ..., xn,yn,...]
        //! @param [in]  n_neigh Number of nearest neighbors of each query point
        //! @param [out] idx     Nearest neighbor list: [n0_0,n0_1,..., n1_0,n1_1,..., ..., nn_0,nn_1,...]
        //! @param [out] dist_sq Nearest distance square: [d0_0,d0_1,..., d1_0,d1_1,..., ..., dn_0,dn_1,...]
        //! @param [in]  eps     Tolerance for nearest neighbor searching. Default is ZERO, which means precision nearest.
        void search_omp(size_type n_pnts, const coord_type* pnts, size_type n_neigh, size_type* idxs, coord_type* dist_sq, coord_type eps = 0)const
        {
            ASSERT(n_neigh <= static_cast<size_type>(pidx_.size()));
            n_neigh = std::min(n_neigh, static_cast<size_type>(pidx_.size()));
            if (n_pnts <= 0 || n_neigh < 0)return;

#pragma omp parallel
            {
                int nt = omp_get_num_threads(); // number of threads
                int it = omp_get_thread_num();  // id of current thread
                int dt = (n_pnts + (size_type)nt - 1) / nt; // task number of each thread

                KDQuery Q{
                    pnts + ND * it * dt,         // q
                    n_neigh,                     // n_neigh
                    0,                           // n_found
                    dist_sq + it * dt * n_neigh, // key
                    idxs + it * dt * n_neigh,    // info
                    (1.0 + eps) * (1.0 + eps),   // max_err
                    0                            // n_visited
                };

                int end = std::min((it + 1) * dt, n_pnts);
                for (int i = dt * it; i < end; ++i) {
                    // compute distance from point to box
                    coord_type box_dist = 0;
                    for (size_t d = 0; d < ND; ++d) {
                        if (Q.q[d] < bnd_box_lo_[d]) {
                            auto delt = bnd_box_lo_[d] - Q.q[d];
                            box_dist += delt * delt;
                        }
                        else if (Q.q[d] > bnd_box_hi_[d]) {
                            auto delt = bnd_box_hi_[d] - Q.q[d];
                            box_dist += delt * delt;
                        }
                    }
                    Q.n_found = 0;

                    // search
                    search_(root_, Q, box_dist);

                    // fill other
                    for (size_type j = Q.n_found; j < Q.n_neigh; ++j) {
                        Q.key[j] = INF;
                        Q.info[j] = -1;
                    }

                    Q.q += ND;
                    Q.n_found = 0;
                    Q.key += n_neigh;
                    Q.info += n_neigh;
                    Q.n_visited = 0;
                }
            }
        }

        //! @brief Is this KDTree empty?
        bool empty()const noexcept { return pidx_.empty(); }

        //! @brief get number of points
        //! @return 
        size_type size()const noexcept { return static_cast<size_type>(pidx_.size()); }

        //! @brief get address of the coordinates array
        const coord_type* data()const { return reinterpret_cast<const coord_type*>(points_); }

        //! @brief get lower of the bounding box
        const point_type& bbox_lo()const noexcept { return bnd_box_lo_; }

        //! @brief get upper of the bounding box
        const point_type& bbox_hi()const noexcept { return bnd_box_hi_; }

        //! @brief get point coordinates
        const point_type& point(size_type i)const { ASSERT(i >= 0 && i < pidx_.size()); return points_[i]; }

        //! @brief get dimension of the KDTree
        constexpr size_t ndim()const noexcept { return NDIM; }

    private:

        struct KDNode
        {
            const point_type* point{ nullptr };
            size_t            cut_dim{ 0 };
            coord_type        cut_val{ 0 };
            coord_type        cd_bnds[2]{ 0 };
            KDNode* child[2]{ nullptr };
        };

        struct KDQuery
        {
            const coord_type* q;
            size_type         n_neigh;
            size_type         n_found;
            coord_type* key;
            size_type* info;
            coord_type        max_err;
            int               n_visited;
        };

        void init_(const coord_type* coords, size_type npts, bool persistent)
        {
            // coordinates
            if (persistent) {
                points_ = reinterpret_cast<const point_type*>(coords);
                points_stg_.clear();
            }
            else {
                points_stg_.resize(npts);
                std::memcpy(points_stg_.data(), coords, sizeof(coord_type) * ND * npts);
                points_ = points_stg_.data();
            }

            // point index
            pidx_.resize(npts);
            for (size_type i = 0; i < npts; ++i)pidx_[i] = i;

            // nodes
            nodes_.clear();
            nodes_.reserve(2 * (size_t)npts); //? 2*npts is sufficient!!!
            root_ = nullptr;

            // bounding box
            bnd_box_lo_.fill(0);
            bnd_box_hi_.fill(0);
        }

        void bounding_box_()noexcept
        {
            if (pidx_.empty())return;

            bnd_box_lo_ = bnd_box_hi_ = points_[0];
            for (size_t i = 1; i < pidx_.size(); ++i) {
                auto& p = points_[i];
                for (size_t j = 0; j < ND; ++j) {
                    auto c = p[j];
                    if (c < bnd_box_lo_[j])bnd_box_lo_[j] = c;
                    else if (c > bnd_box_hi_[j])bnd_box_hi_[j] = c;
                }
            }
        }

        KDNode* build_(size_type* pidx, size_type n, point_type& bnd_box_lo, point_type& bnd_box_hi)
        {
            if (n == 0)return nullptr;
            else if (n == 1) {
                nodes_.emplace_back(KDNode{});
                KDNode* node = &nodes_.back();// nodes_ + nused_; ++nused_; assert(nused_ <= nnode_);
                node->point = points_ + *pidx;
                return node;
            }
            else {
                size_t cd = 0;
                coord_type cv = 0;
                size_type n_lo = 0;

                // find length of longest box side
                coord_type max_length{ 0 };
                for (size_t d = 0; d < ND; ++d) {
                    max_length = std::max(bnd_box_hi[d] - bnd_box_lo[d], max_length);
                }

                // find long side with most spread
                coord_type max_spread = -1;
                coord_type cmin = 0, cmax = 0;
                for (size_t d = 0; d < ND; ++d) {
                    // is it among longest?
                    if ((bnd_box_hi[d] - bnd_box_lo[d]) >= (coord_type{ 1 } - ERR) * max_length) {
                        // compute its spread
                        //coord_type spr = spread_();

                        // compute max and min coords
                        auto xmin = points_[*pidx][d];
                        auto xmax = xmin;
                        for (size_type i = 1; i < n; ++i) {
                            auto c = points_[pidx[i]][d];
                            if (c < xmin)xmin = c;
                            else if (c > xmax)xmax = c;
                        }
                        // total spread is difference
                        coord_type spr = xmax - xmin;

                        // is it max so far?
                        if (spr > max_spread) {
                            max_spread = spr;
                            cd = d;
                            cmin = xmin;
                            cmax = xmax;
                        }
                    }
                }

                // ideal split at midpoint
                auto ideal_cut_val = (bnd_box_lo[cd] + bnd_box_hi[cd]) / 2;

                // find min/max coordinates
                //coord_type cmin, cmax;
                if (max_spread == -1) {
                    cmin = cmax = points_[*pidx][cd];
                    for (size_type i = 1; i < n; ++i) {
                        auto x = points_[pidx[i]][cd];
                        if (x < cmin)cmin = x;
                        else if (x > cmax)cmax = x;
                    }
                }

                // slide to min or max as needed
                if (ideal_cut_val < cmin)
                    cv = cmin;
                else if (ideal_cut_val > cmax)
                    cv = cmax;
                else
                    cv = ideal_cut_val;

                // permute points accordingly
                const ptrdiff_t np = n;
                ptrdiff_t l = 0;
                ptrdiff_t r = np - 1;
                for (;;) {   // partition pa[0..n-1] about cv
                    while (l < np && points_[pidx[l]][cd] < cv) ++l;
                    while (r >= 0 && points_[pidx[r]][cd] >= cv) --r;
                    if (l > r) break;
                    std::swap(pidx[l], pidx[r]);
                    ++l; --r;
                }
                auto br1 = l; // now: pa[0..br1-1] < cv <= pa[br1..n-1]
                r = np - 1;
                for (;;) {							// partition pa[br1..n-1] about cv
                    while (l < np && points_[pidx[l]][cd] <= cv) ++l;
                    while (r >= br1 && points_[pidx[r]][cd] > cv) --r;
                    if (l > r) break;
                    std::swap(pidx[l], pidx[r]);
                    ++l; --r;
                }
                auto br2 = l; // now: pa[br1..br2-1] == cv < pa[br2..n-1]

                //------------------------------------------------------------------
                //	Now:	pa[0....br1-1] <  cut_val
                //			pa[br1..br2-1] == cut_val
                //			pa[br2....n-1] >  cut_val
                //
                //	We can set n_lo to any value in the range [br1..br2] to satisfy
                //	the exit conditions of the procedure.
                //
                //	if ideal_cut_val < min (implying br2 >= 1),
                //			then we select n_lo = 1 (so there is one point on left) and
                //	if ideal_cut_val > max (implying br1 <= n-1),
                //			then we select n_lo = n-1 (so there is one point on right).
                //	Otherwise, we select n_lo as close to n/2 as possible within
                //			[br1..br2].
                //------------------------------------------------------------------
                if (ideal_cut_val < cmin)   n_lo = 1;
                else if (ideal_cut_val > cmax)   n_lo = n - 1;
                else if ((size_type)br1 > n / 2) n_lo = (size_type)br1;
                else if ((size_type)br2 < n / 2) n_lo = (size_type)br2;
                else                             n_lo = n / 2;

                // save bounds for cutting dimension
                auto lv = bnd_box_lo[cd];
                auto hv = bnd_box_hi[cd];

                // modify bounds for left subtree
                bnd_box_hi[cd] = cv;

                // build left subtree: from pidx[0...n_lo-1]
                auto lo = build_(pidx, n_lo, bnd_box_lo, bnd_box_hi);

                // restore bounds
                bnd_box_hi[cd] = hv;

                // modify bounds for right subtree
                bnd_box_lo[cd] = cv;

                // build right subtree
                auto hi = build_(pidx + n_lo, n - n_lo, bnd_box_lo, bnd_box_hi);

                // restore bounds
                bnd_box_lo[cd] = lv;

                // create the splitting node
                nodes_.emplace_back(KDNode{});
                KDNode* node = &nodes_.back();// nodes_ + nused_; ++nused_; assert(nused_ <= nnode_);
                node->cut_dim = cd;
                node->cut_val = cv;
                node->cd_bnds[LO] = lv;
                node->cd_bnds[HI] = hv;
                node->child[LO] = lo;
                node->child[HI] = hi;

                return node;
            }
        }

        void search_(KDNode* node, KDQuery& Q, coord_type box_dist)const
        {
            // search leaf
            if (node->point) {
                ++Q.n_visited;

                // k-th smallest distance so far
                coord_type min_dist = Q.n_found == Q.n_neigh ? Q.key[Q.n_neigh - 1] : INF; //min_dist = ANNkdPointMK->max_key();

                // check points in bucket
                auto pp = node->point->data();// first coord of next data point
                coord_type dist = 0;
                size_t d;
                for (d = 0; d < ND; ++d) {
                    coord_type t = Q.q[d] - pp[d];
                    dist += t * t;
                    if (dist > min_dist)break;
                }

                if (d >= ND /*&& dist != 0*/) {
                    // insert
                    // slide larger values up
                    // +---+---+-...-+--------+
                    // 0   1   2     nfound-1 n-1
                    size_type i = Q.n_found;
                    for (; i > 0; --i) {
                        if (Q.key[i - 1] > dist) {
                            if (i < Q.n_neigh) {
                                Q.key[i] = Q.key[i - 1];
                                Q.info[i] = Q.info[i - 1];
                            }
                        }
                        else break;
                    }
                    if (i < Q.n_neigh) {
                        Q.key[i] = dist;
                        Q.info[i] = static_cast<size_type>(node->point - points_);
                    }
                    if (Q.n_found < Q.n_neigh)++Q.n_found;

                    min_dist = Q.key[Q.n_found - 1];
                }

                return;
            }
            // search split
            else [[likely]] {
                // check dist calc term condition
                // if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

                // distance to cutting plane
                auto cut_diff = Q.q[node->cut_dim] - node->cut_val;
                // left of cutting plane
                if (cut_diff < 0) {
                    // visit closer child first
                    search_(node->child[LO], Q, box_dist);

                    auto box_diff = node->cd_bnds[LO] - Q.q[node->cut_dim];
                    if (box_diff < 0)box_diff = 0;// within bounds - ignore

                    // distance to further box
                    box_dist += cut_diff * cut_diff - box_diff * box_diff;

                    // visit further child if close enough
                    if (box_dist * Q.max_err < (Q.n_found == Q.n_neigh ? Q.key[Q.n_neigh - 1] : INF)) {
                        search_(node->child[HI], Q, box_dist);
                    }
                }
                // right of cutting plane
                else {
                    // visit closer child first
                    search_(node->child[HI], Q, box_dist);

                    auto box_diff = Q.q[node->cut_dim] - node->cd_bnds[HI];
                    if (box_diff < 0)box_diff = 0;// within bounds - ignore

                    // distance to further box
                    box_dist += cut_diff * cut_diff - box_diff * box_diff;

                    // visit further child if close enough
                    if (box_dist * Q.max_err < (Q.n_found == Q.n_neigh ? Q.key[Q.n_neigh - 1] : INF)) {
                        search_(node->child[LO], Q, box_dist);
                    }
                }
            }
        }

        void copy_from_(const KDTree& tree)
        {
            if (&tree == this)return;

            points_ = tree.points_;
            pidx_ = tree.pidx_;
            nodes_ = tree.nodes_;
            bnd_box_lo_ = tree.bnd_box_lo_;
            bnd_box_hi_ = tree.bnd_box_hi_;
            points_stg_ = tree.points_stg_;
            root_ = nullptr;

            // setup pointers
            if (tree.root_) {
                auto n0 = nodes_.data();
                root_ = n0 + (tree.root_ - tree.nodes_.data());

                auto it_p = nodes_.begin();
                auto it_q = tree.nodes_.begin();
                for (; it_p != nodes_.end(); ++it_p, ++it_p) {
                    if (it_q->point)
                        it_p->point = points_ + (it_q->point - tree.points_);
                    else {
                        it_p->child[LO] = n0 + (it_q->child[LO] - tree.nodes_.data());
                        it_p->child[HI] = n0 + (it_q->child[HI] - tree.nodes_.data());
                    }
                }
            }
        }

        void move_(KDTree&& tree)noexcept
        {
            if (&tree == this)return;

            points_ = tree.points_;
            pidx_ = std::move(tree.pidx_);
            nodes_ = std::move(tree.nodes_);
            root_ = tree.root_;
            bnd_box_lo_ = tree.bnd_box_lo_;
            bnd_box_hi_ = tree.bnd_box_hi_;
            points_stg_ = tree.points_stg_;

            tree.points_ = nullptr;
            tree.root_ = nullptr;
            tree.bnd_box_lo_.fill(0);
            tree.bnd_box_hi_.fill(0);
        }

    private:
        const point_type* points_{ nullptr };
        KDNode* root_{ nullptr };
        point_type              bnd_box_lo_{ 0 };
        point_type              bnd_box_hi_{ 0 };
        std::vector<size_type>  pidx_;
        std::vector<KDNode>     nodes_;
        std::vector<point_type> points_stg_;
    };
}
