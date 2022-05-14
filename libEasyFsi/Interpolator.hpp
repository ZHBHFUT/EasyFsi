#pragma once
#include <vector>

namespace EasyLib {

    class Boundary;
    struct Field;

    class Interpolator
    {
    public:
        inline static constexpr const int max_donor = 20;

        Interpolator() = default;

        void clear();

        void set_app_id(int source_app, int target_app);

        //! @brief Add a source boundary to interpolator. All fields of source boundary must be node-centered.
        //! @param bd The source boundary
        //! @return index in the interpolator.
        int add_source_boundary(Boundary& bd);

        //! @brief Add a target boundary to interpolator.
        //! @param bd The target boundary
        //! @return index in the interpolator.
        int add_target_boundary(Boundary& bd);

        //! @brief Compute interpolation coefficients.
        void compute_interp_coeff();

        //! @brief Interpolate incoming DOFs for target boundaries.
        //! @param sources Field of source boundaries. size = source bounds number.
        //! @param targets Field of target boundaries. size = target bounds number.
        //! @note each of \sources field must be node-centered.
        void interp_target_dofs (std::span<Field* const> sources, std::span<Field*> targets);

        //! @brief Interpolate incoming LOADs for source boundaries. 
        //! @param targets Field of target boundaries. size = target bounds number.
        //! @param sources Field of source boundaries. size = source bounds number.
        //! @note each of \sources field must be node-centered.
        void interp_source_loads(std::span<Field* const> targets, std::span<Field*> sources);

        inline int source_app_id()const noexcept { return source_app_; }
        inline int target_app_id()const noexcept { return target_app_; }

    private:
        int source_app_{ -1 };
        int target_app_{ -1 };
        std::vector<Boundary*> source_bounds_;
        std::vector<Boundary*> target_bounds_;
        bool computed_{ false };

        struct Coeffs
        {
            int    src_bd_id{ 0 };
            int    ndonor{ 0 };
            int_l  donor_nodes[max_donor]{ 0 };
            double donor_weights[max_donor]{ 0 };
            double dist_sq{ 0 };
        };
        std::vector<Coeffs> node_coeffs_;
        std::vector<Coeffs> face_coeffs_;
    };
}
