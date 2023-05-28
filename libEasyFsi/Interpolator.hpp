#pragma once
#include <limits>

#include "DynamicMatrix.hpp"

namespace EasyLib {

    class Boundary;
    struct Field;

    enum InterpolationMethod
    {
        //GlobalXPS,  //! Use global IPS/TPS/Spline method.
        LocalXPS,   //! Use local IPS/TPS/Spline method.
        Projection, //! Used projection method, faces should be defined on boundary.
        Mapping,    //! This is not implemented now.
        Automatic   //! Automatic select above method.
    };

    class Interpolator
    {
    public:
        inline static constexpr const int max_donor = 20;
        inline static constexpr const int min_donor = 8;

        Interpolator() = default;
        virtual ~Interpolator() = default;

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
        void compute_interp_coeff(InterpolationMethod method = Automatic, int max_donor_for_xps = max_donor);

        //! @brief Interpolate incoming DOFs from source bounds to target bounds.
        //! @param sources Field of source boundaries. size = source bounds number.
        //! @param targets Field of target boundaries. size = target bounds number.
        //! @note each of \sources field must be node-centered.
        void interp_dofs_s2t(std::span<Field* const> sources, std::span<Field*> targets)const;

        //! @brief Interpolate incoming LOADs from target bounds to source bounds. 
        //! @param targets Field of target boundaries. size = target bounds number.
        //! @param sources Field of source boundaries. size = source bounds number.
        //! @param fill_src_zeros_first  Whether or not to fill the source fields with zeros before computing.
        //! @note
        //!  + each of \sources field must be node-centered.
        void interp_load_t2s(std::span<Field* const> targets, std::span<Field*> sources, bool fill_src_zeros_first = true)const;

        void interp_node_dofs_s2t(int nfield, const double** src_node_dofs, double** des_node_dofs)const;
        void interp_face_dofs_s2t(int nfield, const double** src_node_dofs, double** des_face_dofs)const;

        void interp_node_load_t2s(int nfield, double** src_node_load, const double** des_node_load, bool fill_src_zeros_first = true)const;
        void interp_face_load_t2s(int nfield, double** src_node_load, const double** des_face_load, bool fill_src_zeros_first = true)const;

        inline int source_app_id()const noexcept { return source_app_; }
        inline int target_app_id()const noexcept { return target_app_; }

        void save_coefficients(const char* file)const;
        void read_coefficients(const char* file);

        void write_tecplot(const char* file, int nfield, const char** field_names, const double** src_node_fields, const double** des_node_fields)const;

    private:
        int source_app_{ -1 };
        int target_app_{ -1 };
        std::vector<Boundary*> source_bounds_;
        std::vector<Boundary*> target_bounds_;
        bool computed_{ false };

        struct InterpInfo
        {
            int    src_bd_id{ 0 };
            int    ndonor{ 0 };
            int_l  donor_nodes[max_donor]{ 0 };
            double donor_weights[max_donor]{ 0 };
            size_t donor_beg{ 0 }, donor_end{ 0 };
            double dist_sq{ std::numeric_limits<double>::max() };
        };
        std::vector<InterpInfo> node_info_;
        std::vector<InterpInfo> face_info_;
    };
}
