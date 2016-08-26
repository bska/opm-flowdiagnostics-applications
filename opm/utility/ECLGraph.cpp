/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <opm/utility/ECLGraph.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <exception>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <boost/filesystem.hpp>

#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_nnc_export.h>
#include <ert/util/ert_unique_ptr.hpp>

/// \file
///
/// Implementation of \c ECLGraph interface.

namespace {
    namespace ECL {
        using GridPtr = ::ERT::ert_unique_ptr<ecl_grid_type, ecl_grid_free>;
        using FilePtr = ::ERT::ert_unique_ptr<ecl_file_type, ecl_file_close>;

        /// Internalise on-disk representation of ECLIPSE grid.
        ///
        /// \param[in] grid Name or prefix of on-disk representation of
        ///                 ECLIPSE grid.  If using a prefix, the loader
        ///                 will consider both .EGRID and .GRID versions of
        ///                 the input.
        ///
        /// \return Internalised ERT Grid representation.
        GridPtr loadCase(const boost::filesystem::path& grid);

        /// Internalise on-disk representation of ECL file.
        ///
        /// \param[in] file Name of ECLIPSE output file.
        ///
        /// \return Internalised ERT file contents.
        FilePtr loadFile(const boost::filesystem::path& file);

        /// Retrieve total number of grids managed by model's main grid.
        ///
        /// \param[in] G Main grid obtained from loadCase().
        ///
        /// \return Total number of grids in \p G.  Equal to 1 + total
        /// number of LGRs in model.
        int numGrids(const ecl_grid_type* G);

        /// Access individual grid by numeric index.
        ///
        /// \param[in] G Main grid obtained from loadCase().
        ///
        /// \param[in] gridID Numeric index of requested grid.  Zero for the
        ///    main grid (i.e., \p G itself) or positive for one of the
        ///    LGRs.  Must be strictly less than \code numGrids(G) \endcode.
        ///
        /// \return Pointer to ECL grid corresponding to numeric ID.
        const ecl_grid_type*
        getGrid(const ecl_grid_type* G, const int gridID);

        /// Extract Cartesian dimensions of an ECL grid.
        ///
        /// \param[in] G ERT grid instance corresponding to the model's main
        ///    grid or one of its LGRs.  Typically obtained from function
        ///    getGrid().
        ///
        /// \return Cartesian dimensions of \p G.  Corresponds to number of
        ///    cells in each cardinal direction in 3D depositional space.
        std::array<std::size_t,3>
        cartesianDimensions(const ecl_grid_type* G);

        /// Retrieve global pore-volume vector from INIT source.
        ///
        /// Specialised tool needed to determine the active cells.
        ///
        /// \param[in] G ERT Grid representation.
        ///
        /// \param[in] init ERT representation of INIT source.
        ///
        /// \return Vector of pore-volumes for all global cells of \p G.
        std::vector<double>
        getPVolVector(const ecl_grid_type* G,
                      const ecl_file_type* init,
                      const int            grid_ID = 0);

        /// Extract non-neighbouring connections from ECLIPSE model
        ///
        /// \param[in] G ERT Grid representation corresponding to model's
        ///    main grid obtained directly from loadCase().
        ///
        /// \param[in] init ERT representation of INIT source.
        ///
        /// \return Model's non-neighbouring connections, including those
        ///    between main and local grids.
        std::vector<ecl_nnc_type>
        loadNNC(const ecl_grid_type* G,
                const ecl_file_type* init);

        class CartesianGridData
        {
        public:
            /// Constructor.
            ///
            /// \param[in] G ERT grid structure corresponding either to the
            ///    model's main grid or, if applicable, one of its LGRs.
            ///
            /// \param[in] init Internalised ERT representation of result
            ///    set's INIT file.
            ///
            /// \param[in] gridID Numeric identifier of this grid.  Zero for
            ///    main grid, positive for LGRs.
            CartesianGridData(const ecl_grid_type* G,
                              const ecl_file_type* init,
                              const int            gridID);

            /// Retrieve number of active cells in graph.
            std::size_t numCells() const;

            /// Retrive number of connections in graph.
            std::size_t numConnections() const;

            /// Retrive neighbourship relations between active cells.
            ///
            /// The \c i-th connection is between active cells \code
            /// neighbours()[2*i + 0] \endcode and \code neighbours()[2*i + 1]
            /// \endcode.
            const std::vector<int>& neighbours() const;

            /// Retrive static pore-volume values on active cells only.
            ///
            /// Corresponds to the \c PORV vector in the INIT file, possibly
            /// restricted to those active cells for which the pore-volume is
            /// strictly positive.
            const std::vector<double>& activePoreVolume() const;

            /// Retrieve ID of active cell from global ID.
            int activeCell(const std::size_t globalCell) const;

            /// Retrieve ID of active cell from (I,J,K) index tuple.
            int activeCell(const int i, const int j, const int k) const;

            /// Predicate for whether or not a particular active cell is
            /// further subdivided by an LGR.
            ///
            /// \param[in] cellID Index of particular active cell in this
            ///     grid.
            ///
            /// \return Whether or not cell identified by grid-local active
            ///     index \p cellID is further subdivided by an LGR.
            bool isSubdivided(const int cellID) const;

            /// Retrieve values of result set vector for all global cells in
            /// grid.
            ///
            /// Mostly for implementing connectionData().
            ///
            /// \param[in] src ECLIPSE result set.
            ///
            /// \param[in] vector Name of result set vector.
            ///
            /// \return Numerical values of result set vector, relative to
            /// global cell numbering of this grid.
            std::vector<double>
            cellData(const ecl_file_type* src,
                     const std::string&   vector) const;

            /// Retrieve values of result set vector for all Cartesian
            /// connections in grid.
            ///
            /// \param[in] src ECLIPSE result set.
            ///
            /// \param[in] vector Name prefix of result set vector (e.g.,
            ///     "FLROIL" for oil flux (flow-rate of oil)).
            ///
            /// \return Numerical values of result set vector attributed to
            ///     all of the grid's Cartesian connections.
            std::vector<double>
            connectionData(const ecl_file_type* src,
                           const std::string&   vector) const;

        private:
            /// Facility for deriving Cartesian neighbourship in a grid
            /// (main or LGR) and for mapping result set vectors to grid's
            /// canonical (global) cells.
            class CartesianCells
            {
            public:
                /// Canonical directions of Cartesian neighbours.
                enum class Direction { I, J, K };

                /// Constructor
                ///
                /// \param[in] G ERT Grid representation.
                ///
                /// \param[in] pvol Vector of pore-volumes on all global
                ///                 cells of \p G.  Typically obtained
                ///                 through function getPVolVector().
                CartesianCells(const ecl_grid_type*       G,
                               const std::vector<double>& pvol);

                /// Retrive global cell indices of all active cells in grid.
                std::vector<std::size_t> activeGlobal() const;

                const std::vector<double>& activePoreVolume() const;

                /// Map input vector to all global cells.
                ///
                /// \param[in] x Input vector, defined on the explicitly
                ///              active cells, all global cells or some
                ///              other subset (e.g., all non-neighbouring
                ///              connections).
                ///
                /// \return Input vector mapped to global cells or unchanged
                /// if input is defined on some other subset.
                template <typename T>
                std::vector<T>
                scatterToGlobal(const std::vector<T>& x) const;

                /// Retrieve total number of cells in grid, including
                /// inactive ones.
                ///
                /// Needed to allocate result vectors on global cells.
                std::size_t numGlobalCells() const;

                /// Retrieve active cell ID of particular global cell.
                ///
                /// \param[in] globalCell Index of particular global cell.
                ///
                /// \return Active cell ID of \p globalCell.  Returns
                /// negative one (\code -1 \endcode) if \code globalCell >=
                /// numGlobalCells \endcode or if the global cell is
                /// inactive.
                int getActiveCell(const std::size_t globalCell) const;

                /// Retrieve global cell ID of from (I,J,K) index tuple.
                std::size_t
                getGlobalCell(const int i, const int j, const int k) const;

                /// Retrieve active cell ID of particular global cell's
                /// neighbour in given Cartesian direction.
                ///
                /// \param[in] globalCell Index of particular global cell.
                ///
                /// \param[in] d Cartesian direction in which to look for a
                /// neighbouring cell.
                ///
                /// \return Active cell ID of \p globalCell's neighbour in
                /// direction \d.  Returns negative one (\code -1 \endcode)
                /// if \code globalCell >= numGlobalCells \endcode or if the
                /// global cell is inactive, or if there is no neighbour in
                /// direction \p d (e.g., if purported neighbour would be
                /// outside model).
                int getNeighbour(const std::size_t globalCell,
                                 const Direction   d) const;

                /// Predicate for whether or not a particular active cell is
                /// further subdivided by an LGR.
                bool isSubdivided(const int cellID) const;

            private:
                struct ResultSetMapping {
                    /// Explicit mapping between ACTNUM!=0 cells and global
                    /// cells.
                    struct ID {
                        std::size_t act;
                        std::size_t glob;
                    };

                    /// Number of explicitly active cells (SUM(ACTNUM != 0)).
                    std::size_t num_active;

                    /// Active subset of global cells.
                    std::vector<ID> subset;
                };

                using IndexTuple = std::array<std::size_t,3>;

                /// Size of grid's bounding box (i.e., the number of cells
                /// in each cardinal direction in 3D depositional space).
                const IndexTuple cartesianSize_;

                /// Map cell-based data vectors to grid's global cells.
                ResultSetMapping rsMap_;

                /// Static pore-volumes of all active cells.
                std::vector<double> activePVol_;

                /// Active index of model's global cells.
                std::vector<int> active_ID_;

                /// Whether or not a particular active cell is subdivided.
                std::vector<bool> is_divided_;

                /// Identify those grid cells that are further subdivided by
                /// an LGR.
                ///
                /// Writes to \c is_divided_.
                ///
                /// \param
                void identifySubdividedCells(const ecl_grid_type* G);

                /// Compute linear index of global cell from explicit
                /// (I,J,K) tuple.
                ///
                /// \param[in] ijk Explicit (I,J,K) tuple of global cell.
                ///
                /// \return Linear index (natural ordering) of global cell
                /// (I,J,K).
                std::size_t globIdx(const IndexTuple& ijk) const;

                /// Decompose global (linear) cell index into its (I,J,K)
                /// index tuple.
                ///
                /// \param[in] globalCell Index of particular global cell.
                ///    Must be in the range \code [0 .. numGlobalCells())
                ///    \endcode.
                ///
                /// \return Index triplet of \p globalCell's location within
                /// model.
                IndexTuple ind2sub(const std::size_t globalCell) const;
            };

            /// Collection of global (cell) IDs.
            using GlobalIDColl = std::vector<std::size_t>;

            /// Collection of (global) cell IDs corresponding to the flow
            /// source of each connection.
            using OutCell =
                std::map<CartesianCells::Direction, GlobalIDColl>;

            /// Collection of direction strings to simplify vector name
            /// derivation (replaces chains of if-else)
            using DirectionSuffix =
                std::map<CartesianCells::Direction, std::string>;

            /// Numeric identity of this grid.  Zero for main grid, greater
            /// than zero for LGRs.
            const int gridID_;

            /// Map results from active to global cells.
            CartesianCells cells_;

            /// Known directional suffixes.
            DirectionSuffix suffix_;

            /// Flattened neighbourship relation (array of size \code
            /// 2*numConnections() \endcode).
            std::vector<int> neigh_;

            /// Source cells for each Cartesian connection.
            OutCell outCell_;

            /// Predicate for whether or not a particular result vector is
            /// defined on the grid's cells.
            ///
            /// \param[in] src Result set.
            ///
            /// \param[in] vector Name of result vector.
            ///
            /// \return Whether or not \p vector is defined on model's
            /// cells and part of the result set \p src.
            bool haveCellData(const ecl_file_type* src,
                              const std::string&   vector) const;

            /// Predicate for whether or not a particular result vector is
            /// defined on the grid's Cartesian connections.
            ///
            /// \param[in] src Result set.
            ///
            /// \param[in] vector Prefix of result vector name.
            ///
            /// \return Whether or not all vectors formed by \p vector plus
            /// known directional suffixes are defined on model's cells and
            /// part of the result set \p src.
            bool haveConnData(const ecl_file_type* src,
                              const std::string&   vector) const;

            /// Append directional cell data to global collection of
            /// connection data identified by vector name prefix.
            ///
            /// \param[in] src Result set.
            ///
            /// \param[in] d Cartesian direction.
            ///
            /// \param[in] vector Prefix of result vector name.
            ///
            /// \param[in,out] x Global collection of connection data.  On
            /// input, collection of values corresponding to any previous
            /// directions (preserved), and on output additionally contains
            /// the data corresponding to the Cartesian direction \p d.
            void connectionData(const ecl_file_type*            src,
                                const CartesianCells::Direction d,
                                const std::string&              vector,
                                std::vector<double>&            x) const;

            /// Form complete name of directional result set vector from
            /// prefix and identified direction.
            ///
            /// \param[in] vector Prefix of result vector name.
            ///
            /// \param[in] d Cartesian direction.
            ///
            /// \return \code vector + suffix_[d] \endcode.
            std::string
            vectorName(const std::string&              vector,
                       const CartesianCells::Direction d) const;

            /// Derive neighbourship relations on active cells in particular
            /// Cartesian directions.
            ///
            /// Writes to \c neigh_ and \c outCell_.
            ///
            /// \param[in] gcells Collection of global (relative to \c
            ///    gridID_) cells that should be considered active (strictly
            ///    positive pore-volume and not deactivated through
            ///    ACTNUM=0).
            ///
            /// \param[in] init Internalised
            void deriveNeighbours(const std::vector<std::size_t>& gcells,
                                  const ecl_file_type*            init,
                                  const CartesianCells::Direction d);
        };
    } // namespace ECL
} // Anonymous namespace

// ======================================================================

int ECL::numGrids(const ecl_grid_type* G)
{
    return 1 + ecl_grid_get_num_lgr(G); // Main + #LGRs.
}

const ecl_grid_type*
ECL::getGrid(const ecl_grid_type* G, const int gridID)
{
    assert ((gridID >= 0) && "Grid ID must be non-negative");

    if (gridID == 0) {
        return G;
    }

    return ecl_grid_iget_lgr(G, gridID - 1);
}

std::vector<double>
ECL::getPVolVector(const ecl_grid_type* G,
                   const ecl_file_type* init,
                   const int            gridID)
{
    auto make_szt = [](const int i)
    {
        return static_cast<std::vector<double>::size_type>(i);
    };

    const auto nglob = make_szt(ecl_grid_get_global_size(G));

    auto pvol = std::vector<double>(nglob, 1.0);

    if (ecl_file_has_kw(init, "PORV")) {
        auto porv =
            ecl_file_iget_named_kw(init, "PORV", gridID);

        assert ((make_szt(ecl_kw_get_size(porv)) == nglob)
                && "Pore-volume must be provided for all global cells");

        ecl_kw_get_data_as_double(porv, pvol.data());
    }

    return pvol;
}

ECL::GridPtr
ECL::loadCase(const boost::filesystem::path& grid)
{
    auto G = GridPtr{
        ecl_grid_load_case(grid.generic_string().c_str())
    };

    if (! G) {
        std::ostringstream os;

        os << "Failed to load ECL Grid from '"
           << grid.generic_string() << '\'';

        throw std::invalid_argument(os.str());
    }

    return G;
}

ECL::FilePtr
ECL::loadFile(const boost::filesystem::path& file)
{
    // Read-only, keep open between requests
    const auto open_flags = 0;

    auto F = FilePtr{
        ecl_file_open(file.generic_string().c_str(), open_flags)
    };

    if (! F) {
        std::ostringstream os;

        os << "Failed to load ECL File object from '"
           << file.generic_string() << '\'';

        throw std::invalid_argument(os.str());
    }

    return F;
}

std::array<std::size_t,3>
ECL::cartesianDimensions(const ecl_grid_type* G)
{
    auto make_szt = [](const int i)
    {
        return static_cast<std::size_t>(i);
    };

    return { { make_szt(ecl_grid_get_nx(G)) ,
               make_szt(ecl_grid_get_ny(G)) ,
               make_szt(ecl_grid_get_nz(G)) } };
}

std::vector<ecl_nnc_type>
ECL::loadNNC(const ecl_grid_type* G,
             const ecl_file_type* init)
{
    auto make_szt = [](const int n)
    {
        return static_cast<std::vector<ecl_nnc_type>::size_type>(n);
    };

    auto nncData = std::vector<ecl_nnc_type>{};

    const auto numNNC = make_szt(ecl_nnc_export_get_size(G));

    if (numNNC > 0) {
        nncData.resize(numNNC);

        ecl_nnc_export(G, init, nncData.data());
    }

    return nncData;
}

// ======================================================================

ECL::CartesianGridData::
CartesianCells::CartesianCells(const ecl_grid_type*       G,
                               const std::vector<double>& pvol)
    : cartesianSize_(::ECL::cartesianDimensions(G))
{
    if (pvol.size() != static_cast<decltype(pvol.size())>
        (this->cartesianSize_[0] *
         this->cartesianSize_[1] *
         this->cartesianSize_[2]))
    {
        throw std::invalid_argument("Grid must have PORV for all cells");
    }

    auto make_szt = [](const int i)
    {
        return static_cast<std::size_t>(i);
    };

    using ID = ResultSetMapping::ID;

    this->rsMap_.num_active = make_szt(ecl_grid_get_nactive(G));

    this->rsMap_.subset.clear();
    this->rsMap_.subset.reserve(this->rsMap_.num_active);

    for (decltype(ecl_grid_get_nactive(G))
             act = 0, nact = ecl_grid_get_nactive(G);
         act < nact; ++act)
    {
        const auto glob =
            make_szt(ecl_grid_get_global_index1A(G, act));

        if (pvol[glob] > 0.0) {
            this->rsMap_.subset.push_back(ID{ make_szt(act), glob });
        }
    }

    {
        std::vector<int>(pvol.size(), -1).swap(this->active_ID_);

        this->activePVol_.clear();
        this->activePVol_.reserve(this->rsMap_.subset.size());

        this->is_divided_.clear();
        this->is_divided_.reserve(this->rsMap_.subset.size());

        auto active = 0;

        for (const auto& cell : this->rsMap_.subset) {
            this->active_ID_[cell.glob] = active++;
            this->activePVol_.push_back(pvol[cell.glob]);

            const auto ert_active = static_cast<int>(cell.act);
            const auto is_divided =
                nullptr != ecl_grid_get_cell_lgr1A(G, ert_active);

            this->is_divided_.push_back(is_divided);
        }
    }
}

std::vector<std::size_t>
ECL::CartesianGridData::CartesianCells::activeGlobal() const
{
    auto active = std::vector<std::size_t>{};
    active.reserve(this->rsMap_.subset.size());

    for (const auto& id : this->rsMap_.subset) {
        active.push_back(id.glob);
    }

    return active;
}

const std::vector<double>&
ECL::CartesianGridData::CartesianCells::activePoreVolume() const
{
    return this->activePVol_;
}

template <typename T>
std::vector<T>
ECL::CartesianGridData::
CartesianCells::scatterToGlobal(const std::vector<T>& x) const
{
    // Assume that input vector 'x' is either defined on explicit notion of
    // active cells (ACTNUM != 0) or on all global cells or some other
    // contiguous index set (e.g., the NNCs).

    const auto num_explicit_active =
        static_cast<decltype(x.size())>(this->rsMap_.num_active);

    if (x.size() != num_explicit_active) {
        // Input not defined on explictly active cells.  Let caller deal
        // with it.  This typically corresponds to the set of global cells
        // or the list of NNCs.
        return x;
    }

    auto y = std::vector<T>(this->numGlobalCells());

    for (const auto& i : this->rsMap_.subset) {
        y[i.glob] = x[i.act];
    }

    return y;
}

std::size_t
ECL::CartesianGridData::CartesianCells::numGlobalCells() const
{
    return this->active_ID_.size();
}

int
ECL::CartesianGridData::
CartesianCells::getActiveCell(const std::size_t globalCell) const
{
    if (globalCell >= numGlobalCells()) { return -1; }

    return this->active_ID_[globalCell];
}

std::size_t
ECL::CartesianGridData::
CartesianCells::getGlobalCell(const int i, const int j, const int k) const
{
    const auto ijk = IndexTuple {
        static_cast<std::size_t>(i),
        static_cast<std::size_t>(j),
        static_cast<std::size_t>(k),
    };

    return this->globIdx(ijk);
}

int
ECL::CartesianGridData::
CartesianCells::getNeighbour(const std::size_t globalCell,
                             const Direction   d) const
{
    if (globalCell >= numGlobalCells()) { return -1; }

    auto ijk = ind2sub(globalCell);

    if      (d == Direction::I) { ijk[0] += 1; }
    else if (d == Direction::J) { ijk[1] += 1; }
    else if (d == Direction::K) { ijk[2] += 1; }
    else {
        return -1;
    }

    const auto globNeigh = globIdx(ijk);

    if (globNeigh >= numGlobalCells()) { return -1; }

    return this->active_ID_[globNeigh];
}

bool
ECL::CartesianGridData::CartesianCells::isSubdivided(const int cellID) const
{
    const auto ix =
        static_cast<decltype(this->is_divided_.size())>(cellID);

    assert ((cellID >= 0) && (ix < this->is_divided_.size()));

    return this->is_divided_[ix];
}

std::size_t
ECL::CartesianGridData::
CartesianCells::globIdx(const IndexTuple& ijk) const
{
    const auto& dim = this->cartesianSize_;

    for (auto d = 0*dim.size(), nd = dim.size(); d < nd; ++d) {
        if (ijk[d] >= dim[d]) { return -1; }
    }

    return ijk[0] + dim[0]*(ijk[1] + dim[1]*ijk[2]);
}

ECL::CartesianGridData::CartesianCells::IndexTuple
ECL::CartesianGridData::
CartesianCells::ind2sub(const std::size_t globalCell) const
{
    assert (globalCell < numGlobalCells());

    auto ijk = IndexTuple{};
    auto g   = globalCell;

    const auto& dim = this->cartesianSize_;

    ijk[0] = g % dim[0];  g /= dim[0];
    ijk[1] = g % dim[1];
    ijk[2] = g / dim[1];  assert (ijk[2] < dim[2]);

    assert (globIdx(ijk) == globalCell);

    return ijk;
}

// ======================================================================

ECL::CartesianGridData::CartesianGridData(const ecl_grid_type* G,
                                          const ecl_file_type* init,
                                          const int            gridID)
    : gridID_(gridID)
    , cells_ (G, ::ECL::getPVolVector(G, init, gridID_))
{
    {
        using VT = DirectionSuffix::value_type;

        suffix_.insert(VT(CartesianCells::Direction::I, "I+"));
        suffix_.insert(VT(CartesianCells::Direction::J, "J+"));
        suffix_.insert(VT(CartesianCells::Direction::K, "K+"));
    }

    const auto gcells = this->cells_.activeGlobal();

    // Too large, but this is a quick estimate.
    this->neigh_.reserve(3 * (2 * this->numCells()));

    for (const auto d : { CartesianCells::Direction::I ,
                          CartesianCells::Direction::J ,
                          CartesianCells::Direction::K })
    {
        this->deriveNeighbours(gcells, init, d);
    }
}

std::size_t
ECL::CartesianGridData::numCells() const
{
    return this->activePoreVolume().size();
}

std::size_t
ECL::CartesianGridData::numConnections() const
{
    return this->neigh_.size() / 2;
}

const std::vector<int>&
ECL::CartesianGridData::neighbours() const
{
    return this->neigh_;
}

const std::vector<double>&
ECL::CartesianGridData::activePoreVolume() const
{
    return this->cells_.activePoreVolume();
}

int
ECL::CartesianGridData::activeCell(const std::size_t globalCell) const
{
    return this->cells_.getActiveCell(globalCell);
}

int
ECL::CartesianGridData::activeCell(const int i,
                                   const int j,
                                   const int k) const
{
    return this->activeCell(this->cells_.getGlobalCell(i, j, k));
}

bool
ECL::CartesianGridData::isSubdivided(const int cellID) const
{
    return this->cells_.isSubdivided(cellID);
}

std::vector<double>
ECL::CartesianGridData::cellData(const ecl_file_type* src,
                                 const std::string&   vector) const
{
    if (! this->haveCellData(src, vector)) {
        return {};
    }

    const auto v =
        ecl_file_iget_named_kw(src, vector.c_str(), this->gridID_);

    auto x = std::vector<double>(ecl_kw_get_size(v));

    ecl_kw_get_data_as_double(v, x.data());

    return this->cells_.scatterToGlobal(x);
}

bool
ECL::CartesianGridData::haveCellData(const ecl_file_type* src,
                                     const std::string&   vector) const
{
    // Recall: get_num_named_kw() is block aware (uses src->active_map).

    return ecl_file_get_num_named_kw(src, vector.c_str()) > this->gridID_;
}

bool
ECL::CartesianGridData::haveConnData(const ecl_file_type* src,
                                     const std::string&   vector) const
{
    auto have_data = true;

    for (const auto& d : { CartesianCells::Direction::I ,
                           CartesianCells::Direction::J ,
                           CartesianCells::Direction::K })
    {
        const auto vname = this->vectorName(vector, d);
        have_data = this->haveCellData(src, vname);

        if (! have_data) { break; }
    }

    return have_data;
}

std::vector<double>
ECL::CartesianGridData::connectionData(const ecl_file_type* src,
                                       const std::string&   vector) const
{
    if (! this->haveConnData(src, vector)) {
        return {};
    }

    auto x = std::vector<double>{};  x.reserve(this->numConnections());

    for (const auto& d : { CartesianCells::Direction::I ,
                           CartesianCells::Direction::J ,
                           CartesianCells::Direction::K })
    {
        this->connectionData(src, d, vector, x);
    }

    return x;
}

void
ECL::CartesianGridData::
connectionData(const ecl_file_type*            src,
               const CartesianCells::Direction d,
               const std::string&              vector,
               std::vector<double>&            x) const
{
    const auto v = this->cellData(src, this->vectorName(vector, d));

    const auto& cells = this->outCell_.find(d);

    assert ((cells != this->outCell_.end()) &&
            "Direction must be I, J, or K");

    for (const auto& cell : cells->second) {
        x.push_back(v[cell]);
    }
}

std::string
ECL::CartesianGridData::
vectorName(const std::string&              vector,
           const CartesianCells::Direction d) const
{
    const auto i = this->suffix_.find(d);

    assert ((i != this->suffix_.end()) &&
            "Direction must be I, J, or K");

    return vector + i->second;
}

void
ECL::CartesianGridData::
deriveNeighbours(const std::vector<std::size_t>& gcells,
                 const ecl_file_type*            init,
                 const CartesianCells::Direction d)
{
    auto tran = std::string{"TRAN"};

    switch (d) {
    case CartesianCells::Direction::I:
        tran += 'X';
        break;

    case CartesianCells::Direction::J:
        tran += 'Y';
        break;

    case CartesianCells::Direction::K:
        tran += 'Z';
        break;

    default:
        throw std::invalid_argument("Input direction must be (I,J,K)");
    }

    const auto& T = this->haveCellData(init, tran)
        ? this->cellData(init, tran)
        : std::vector<double>(this->cells_.numGlobalCells(), 1.0);

    auto& ocell = this->outCell_[d];
    ocell.reserve(gcells.size());

    for (const auto& globID : gcells) {
        const auto c1 = this->cells_.getActiveCell(globID);

        assert ((c1 >= 0) && "Internal error in active cell derivation");

        if (this->cells_.isSubdivided(c1)) {
            // Don't form connections to subdivided cells.  We care only
            // about the final refinement level (i.e., the most nested LGR
            // object) and the connections are handled by the NNC code.
            continue;
        }

        if (T[globID] > 0.0) {
            const auto other = this->cells_.getNeighbour(globID, d);

            if ((other >= 0) && ! this->cells_.isSubdivided(other)) {
                assert (c1 != other);

                this->neigh_.push_back(c1);
                this->neigh_.push_back(other);

                ocell.push_back(globID);
            }
        }
    }
}

// =====================================================================

/// Implementation of ECLGraph interface.
class Opm::ECLGraph::Impl
{
public:
    /// Constructor
    ///
    /// \param[in] grid Name or prefix of ECL grid (i.e., .GRID or
    ///                 .EGRID) file.
    ///
    /// \param[in] init Name of ECL INIT file corresponding to \p grid
    ///                 input.  Assumed to provide at least a complete set
    ///                 of pore-volume values (i.e., for all global cells
    ///                 defined in the \p grid).
    ///
    ///                 If available in the INIT file, the constructor will
    ///                 also leverage the transmissibility data when
    ///                 constructing the active cell neighbourship table.
    Impl(const Path& grid, const Path& init);

    /// Assign source object for phase flux calculation.
    ///
    /// \param[in] src Name of ECL restart file, possibly unified, from
    ///                which next set of phase fluxes should be retrieved.
    void assignDataSource(const Path& src);

    /// Retrieve active cell ID from (I,J,K) tuple in particular grid.
    ///
    /// \param[in] gridID Identity of specific grid to which to relate the
    ///     (I,J,K) tuple.  Zero for main grid and positive indices for any
    ///     LGRs.  The (I,J,K) indices must be within the ranges implied by
    ///     the specific grid.
    ///
    /// \param[in] ijk Cartesian index tuple of particular cell.
    ///
    /// \return Active ID (relative to linear, global numbering) of cell \p
    ///     ijk from specified grid.  Negative one (-1) if (I,J,K) outside
    ///     valid range or if the specific cell identified by \p ijk and \p
    ///     gridID is not actually active.
    int activeCell(const int                gridID,
                   const std::array<int,3>& ijk) const;

    /// Retrieve number of active cells in graph.
    std::size_t numCells() const;

    /// Retrive number of connections in graph.
    std::size_t numConnections() const;

    /// Retrive neighbourship relations between active cells.
    ///
    /// The \c i-th connection is between active cells \code
    /// neighbours()[2*i + 0] \endcode and \code neighbours()[2*i + 1]
    /// \endcode.
    std::vector<int> neighbours() const;

    /// Retrive static pore-volume values on active cells only.
    ///
    /// Corresponds to the \c PORV vector in the INIT file, possibly
    /// restricted to those active cells for which the pore-volume is
    /// strictly positive.
    std::vector<double> activePoreVolume() const;

    /// Retrive phase flux on all connections defined by \code neighbours()
    /// \endcode.
    ///
    /// Non-"const" because this potentially loads new data from the backing
    /// store into internal cache data structures.
    ///
    /// \param[in] phase Canonical phase for which to retrive flux.
    ///
    /// \param[in] rptstep Selected temporal vector.  Report-step ID.
    ///
    /// \return Flux values corresponding to selected phase and report step.
    /// Empty if unavailable in the result set (e.g., by querying the gas
    /// flux in an oil/water system or if the specified \p occurrence is not
    /// reported due to report frequencies or no flux values are output at
    /// all).
    std::vector<double>
    flux(const BlackoilPhases::PhaseIndex phase,
         const int                        rptstep) const;

private:
    /// Collection of non-Cartesian neighbourship relations attributed to a
    /// particular ECL keyword set (i.e., one of NNC{1,2}, NNC{G,L}, NNCLL).
    class NonNeighKeywordIndexSet
    {
    public:
        /// Establish mapping between particular non-Cartesian neighbourship
        /// relation and particular entry within a grid's keyword data.
        struct Map {
            /// Non-Cartesian neighbourship relation.
            std::size_t neighIdx;

            /// Index into grid's keyword data.
            std::size_t kwIdx;
        };

        using MapCollection = std::vector<Map>;

        /// Record a mapping in a particular grid.
        ///
        /// \param[in] grid Particular model grid for which to record a
        ///     mapping.
        ///
        /// \param[in] entry Individual index map.
        void add(const int grid, Map&& entry);

        /// Retrieve collection of index maps for particular grid.
        ///
        /// \param[in] grid Specific model grid ID.  Must be non-negative
        ///    and strictly less than the total number of grids in the
        ///    model.
        ///
        /// \return Collection of index maps attributable to \p grid.  Empty
        ///    if no such collection exists.
        const MapCollection& getGridCollection(const int grid) const;

    private:
        using KWEntries = std::map<int, MapCollection>;

        /// Collection of all index maps attributable to all grids for this
        /// set of ECL keywords.
        KWEntries subset_;

        /// Return value for the case of no existing collection.
        MapCollection empty_;
    };

    /// Collection of all non-Cartesian neighbourship relations classified
    /// according to connection type.
    class NNC
    {
    public:
        /// Classification of non-Cartesian neighbourship relations.
        enum class Category {
            /// Traditional non-neighbouring connections entirely internal
            /// to a grid.  Typically due to faults or fully unstructured
            /// grid descriptions.  Keywords NNC{1,2}, TRANNNC, and FLR*N+.
            /// Positive from NNC1 to NNC2.
            Normal,

            /// Connections between main grid and LGRs.  Keywords NNCG,
            /// NNCL, TRANGL, and FLR*L+.  Positive from NNCG to NNCL.
            GlobalToLocal,

            /// Connections between LGRs.  Either due to two LGRs being
            /// neighbouring entities in physical space or one LGR being
            /// nested within another.  Keywords NNA{1,2}, TRANLL, and
            /// FLR*A+.  Positive from NNA1 to NNA2.
            Amalgamated,
        };

        /// Map a collection of non-Cartesian neighbourship relations to a
        /// specific flux vector identifier.
        struct FluxRelation {
            /// Flux vector identifier.  Should be one of "N+" for Normal
            /// connections, "L+" for GlobalToLocal connections, and "A+"
            /// for Amalgamated connections.
            std::string fluxID;

            /// Collection of non-Cartesian neighbourship relations.
            NonNeighKeywordIndexSet indexSet;
        };

        /// Constructor.
        NNC();

        /// Potentially record a new non-Cartesian connection.
        ///
        /// Will classify the connection according to the grids involved and
        /// actually record the connection if both cells are active and
        /// neither are subdivided.
        ///
        /// \param[in] grids Collection of all active grids in model.
        ///
        /// \param[in] offset Start index into global linear number for all
        ///    active grids.
        ///
        /// \param[in] nnc Non-neighbouring connection from result set.
        void add(const std::vector<ECL::CartesianGridData>& grids,
                 const std::vector<std::size_t>&            offset,
                 const ecl_nnc_type&                        nnc);

        std::vector<Category> allCategories() const;

        /// Retrieve total number of active non-neighbouring connections.
        std::size_t numConnections() const;

        /// Access all active non-neighbouring connections.
        const std::vector<int>& getNeighbours() const;

        /// Retrieve all non-neighbouring connections of a particular
        /// category (i.e., pertaining to a particular set of keywords).
        ///
        /// \param[in] type Category of non-neighbouring connections.
        ///
        /// \return All non-neighbouring connections of category \p type.
        const FluxRelation& getRelations(const Category& type) const;

    private:
        using KeywordIndexMap = std::map<Category, FluxRelation>;

        /// Active non-Cartesian (non-neighbouring) connections.  Cell IDs
        /// in linear numbering of all model's active cells.
        std::vector<int> neigh_;

        /// Collection of
        KeywordIndexMap keywords_;

        /// Factory for FluxRelations.
        ///
        /// Simplifies implementation of ctor.
        ///
        /// \param[in] cat Class
        FluxRelation makeRelation(const Category cat) const;

        /// Identify connection category from connection's grids.
        ///
        ///   - Normal connection if both grid IDs equal
        ///   - GlobalToLocal if one grid is the main model grid and the
        ///     other is an LGR.
        ///   - Amalgamated if both grids are LGRs.
        ///
        /// \param[in] grid1 Numeric identity of connection's first grid.
        ///     Zero if main grid, positive if LGR.
        ///
        /// \param[in] grid2 Numeric identity of connection's second grid.
        ///     Zero if main grid, positive if LGR.
        ///
        /// \return Category of connection between \p grid1 and \p grid2.
        Category classifyConnection(const int grid1, const int grid2) const;

        /// Check if cell is viable connection endpoint in grid.
        ///
        /// A cell is a viable connection endpoint if it is active within a
        /// particular grid and not further subdivided by an LGR.
        ///
        /// \param[in] grids Collection of all active grids in model.
        ///
        /// \param[in] gridID Numeric identity of connection grid.
        ///     Zero if main grid, positive if LGR.
        ///
        /// \param[in] cellID Global ID (relative to \p grid) of candidate
        ///     connection endpoint.
        ///
        /// \return Whether or not \p cellID is a viable connection endpoint
        ///     within \p gridID.
        bool isViable(const std::vector<ECL::CartesianGridData>& grids,
                      const int         gridID,
                      const std::size_t cellID) const;

        /// Check if connection is viable
        ///
        /// A candidate non-Cartesian connection is viable if both of its
        /// endpoints satisfy the viability criterion.
        ///
        /// \param[in] nnc Candidate non-Cartesian connection.
        ///
        /// \return Whether or not both candidate endpoints satisfy the
        /// viability criterion.
        bool isViable(const std::vector<ECL::CartesianGridData>& grids,
                      const ecl_nnc_type& nnc) const;
    };

    /// Collection of model's non-neighbouring connections--be they within a
    /// grid or between grids.
    NNC nnc_;

    /// Collection of model's grids (main + LGRs).
    std::vector<ECL::CartesianGridData> grid_;

    /// Map each grid's active cellIDs to global numbering (in the index
    /// range \code [0 .. numCells()) \endcode).
    std::vector<std::size_t> activeOffset_;

    /// Current result set.
    ECL::FilePtr src_;

    /// Extract explicit non-neighbouring connections from ECL output.
    ///
    /// Writes to \c neigh_ and \c nncID_.
    ///
    /// \param[in] G ERT Grid representation.
    ///
    /// \param[in] init ERT representation of INIT source.
    ///
    /// \param[in] coll Backing data for neighbourship extraction.
    void defineNNCs(const ecl_grid_type* G,
                    const ecl_file_type* init);

    /// Compute ECL vector basename for particular phase flux.
    ///
    /// \param[in] phase Canonical phase for which to derive ECL vector
    /// basename.
    ///
    /// \return Basename for ECl vector corresponding to particular phase
    /// flux.
    std::string
    flowVector(const BlackoilPhases::PhaseIndex phase) const;

    /// Extract flux values corresponding to particular result set vector
    /// for all identified non-neighbouring connections.
    ///
    /// \param[in] vector Result set vector prefix.  Typically computed by
    ///    method flowVector().
    ///
    /// \param[in,out] flux Numerical values of result set vector.  On
    ///    input, contains all values corresponding to all fully Cartesian
    ///    connections across all active grids.  On output additionally
    ///    contains those values that correspond to the non-neighbouring
    ///    connections (appended onto \p flux).
    void fluxNNC(const std::string&   vector,
                 std::vector<double>& flux) const;
};

// ======================================================================

void
Opm::ECLGraph::Impl::NonNeighKeywordIndexSet::
add(const int grid, Map&& entry)
{
    this->subset_[grid].push_back(std::move(entry));
}

const
Opm::ECLGraph::Impl::NonNeighKeywordIndexSet::MapCollection&
Opm::ECLGraph::Impl::NonNeighKeywordIndexSet::
getGridCollection(const int grid) const
{
    auto coll = this->subset_.find(grid);

    if (coll == this->subset_.end()) {
        // No NNCs of this category for this grid.  Return empty.
        return this->empty_;
    }

    return coll->second;
}

// ======================================================================

Opm::ECLGraph::Impl::NNC::NNC()
{
    using VT = KeywordIndexMap::value_type;

    for (const auto& cat : this->allCategories()) {
        this->keywords_.insert(VT(cat, this->makeRelation(cat)));
    }
}

std::vector<Opm::ECLGraph::Impl::NNC::Category>
Opm::ECLGraph::Impl::NNC::allCategories() const
{
    return { Category::Normal        ,
             Category::GlobalToLocal ,
             Category::Amalgamated   };
}

void
Opm::ECLGraph::Impl::
NNC::add(const std::vector<ECL::CartesianGridData>& grid,
         const std::vector<std::size_t>&            offset,
         const ecl_nnc_type&                        nnc)
{
    if (! this->isViable(grid, nnc)) {
        // At least one endpoint unviable.  Don't record connection.
        return;
    }

    const auto neighIdx = this->numConnections();

    {
        const auto c = grid[nnc.grid_nr1].activeCell(nnc.global_index1);
        const auto o = static_cast<int>(offset[nnc.grid_nr1]);

        this->neigh_.push_back(o + c);
    }

    {
        const auto c = grid[nnc.grid_nr2].activeCell(nnc.global_index2);
        const auto o = static_cast<int>(offset[nnc.grid_nr2]);

        this->neigh_.push_back(o + c);
    }

    const auto cat = this->classifyConnection(nnc.grid_nr1, nnc.grid_nr2);

    auto entry = NonNeighKeywordIndexSet::Map {
        neighIdx,
        static_cast<std::size_t>(nnc.input_index)
    };

    this->keywords_[cat].indexSet.add(nnc.grid_nr2, std::move(entry));
}

std::size_t
Opm::ECLGraph::Impl::NNC::numConnections() const
{
    assert (this->neigh_.size() % 2 == 0);

    return this->neigh_.size() / 2;
}

const std::vector<int>&
Opm::ECLGraph::Impl::NNC::getNeighbours() const
{
    return this->neigh_;
}

const Opm::ECLGraph::Impl::NNC::FluxRelation&
Opm::ECLGraph::Impl::NNC::getRelations(const Category& cat) const
{
    auto r = this->keywords_.find(cat);

    assert ((r != this->keywords_.end()) &&
            "Input category must be Normal, "
            "GlobalToLocal or Amalgamated");

    return r->second;
}

Opm::ECLGraph::Impl::NNC::FluxRelation
Opm::ECLGraph::Impl::NNC::makeRelation(const Category cat) const
{
    switch (cat) {
    case Category::Normal:
        return { "N+", {} };

    case Category::GlobalToLocal:
        return { "L+", {} };

    case Category::Amalgamated:
        return { "A+", {} };
    }

    throw std::invalid_argument("Category must be Normal, "
                                "GlobalToLocal, or Amalgamated");
}

Opm::ECLGraph::Impl::NNC::Category
Opm::ECLGraph::Impl::NNC::
classifyConnection(const int grid1, const int grid2) const
{
    if (grid1 == grid2) {
        return Category::Normal;
    }

    if (grid1 == 0) {           // Main grid
        return Category::GlobalToLocal;
    }

    return Category::Amalgamated;
}

bool
Opm::ECLGraph::Impl::NNC::
isViable(const std::vector<ECL::CartesianGridData>& grids,
         const int                                  gridID,
         const std::size_t                          cellID) const
{
    using GridIndex = decltype(grids.size());
    const auto gIdx = static_cast<GridIndex>(gridID);

    if (gIdx >= grids.size()) {
        return false;
    }

    const auto& G     = grids[gIdx];
    const auto  acell = G.activeCell(cellID);

    return (acell >= 0) && (! G.isSubdivided(acell));
}

bool
Opm::ECLGraph::Impl::NNC::
isViable(const std::vector<ECL::CartesianGridData>& grids,
         const ecl_nnc_type&                        nnc) const
{
    return this->isViable(grids, nnc.grid_nr1, nnc.global_index1)
        && this->isViable(grids, nnc.grid_nr2, nnc.global_index2);
}

// ======================================================================

Opm::ECLGraph::Impl::Impl(const Path& grid, const Path& init)
{
    const auto G = ECL::loadCase(grid);
    auto       I = ECL::loadFile(init);

    const auto numGrids = ECL::numGrids(G.get());

    this->grid_.reserve(numGrids);
    this->activeOffset_.reserve(numGrids + 1);
    this->activeOffset_.push_back(0);

    for (auto gridID = 0*numGrids; gridID < numGrids; ++gridID)
    {
        this->grid_.emplace_back(ECL::getGrid(G.get(), gridID),
                                 I.get(), gridID);

        this->activeOffset_.push_back(this->activeOffset_.back() +
                                      this->grid_.back().numCells());
    }

    this->defineNNCs(G.get(), I.get());
}

void
Opm::ECLGraph::Impl::assignDataSource(const Path& src)
{
    this->src_ = ECL::loadFile(src);
}

int
Opm::ECLGraph::Impl::
activeCell(const int                gridID,
           const std::array<int,3>& ijk) const
{
    const auto gIdx =
        static_cast<decltype(this->grid_.size())>(gridID);

    if (gIdx >= this->grid_.size()) {
        return -1;
    }

    const auto& grid = this->grid_[gIdx];

    const auto active = grid.activeCell(ijk[0], ijk[1], ijk[2]);

    if ((active < 0) || grid.isSubdivided(active)) {
        return -1;
    }

    const auto off = static_cast<int>(this->activeOffset_[gIdx]);

    return off + active;
}

std::size_t
Opm::ECLGraph::Impl::numCells() const
{
    return this->activeOffset_.back();
}

std::size_t
Opm::ECLGraph::Impl::numConnections() const
{
    auto nconn = std::size_t{0};

    for (const auto& G : this->grid_) {
        nconn += G.numConnections();
    }

    return nconn + this->nnc_.numConnections();
}

std::vector<int>
Opm::ECLGraph::Impl::neighbours() const
{
    auto N = std::vector<int>{};

    N.reserve(2 * (this->numConnections() +
                   this->nnc_.numConnections()));

    {
        auto off = this->activeOffset_.begin();

        for (const auto& G : this->grid_) {
            const auto add = static_cast<int>(*off);

            for (const auto& cell : G.neighbours()) {
                N.push_back(cell + add);
            }

            ++off;
        }
    }

    {
        const auto& nnc = this->nnc_.getNeighbours();

        N.insert(N.end(), nnc.begin(), nnc.end());
    }

    return N;
}

std::vector<double>
Opm::ECLGraph::Impl::activePoreVolume() const
{
    auto pvol = std::vector<double>{};
    pvol.reserve(this->numCells());

    for (const auto& G : this->grid_) {
        const auto& pv = G.activePoreVolume();

        pvol.insert(pvol.end(), pv.begin(), pv.end());
    }

    return pvol;
}

std::vector<double>
Opm::ECLGraph::Impl::
flux(const BlackoilPhases::PhaseIndex phase,
     const int                        rptstep) const
{
    if (! ecl_file_has_report_step(this->src_.get(), rptstep)) {
        return {};
    }

    const auto vector = this->flowVector(phase);

    auto v = std::vector<double>{};

    v.reserve(this->numConnections());

    ecl_file_select_rstblock_report_step(this->src_.get(), rptstep);

    for (const auto& G : this->grid_) {
        const auto& q = G.connectionData(this->src_.get(), vector);

        v.insert(v.end(), q.begin(), q.end());
    }

    this->fluxNNC(vector, v);

    return v;
}

void
Opm::ECLGraph::Impl::defineNNCs(const ecl_grid_type* G,
                                const ecl_file_type* init)
{
    for (const auto& nnc : ECL::loadNNC(G, init)) {
        this->nnc_.add(this->grid_, this->activeOffset_, nnc);
    }
}

void
Opm::ECLGraph::Impl::fluxNNC(const std::string&   vector,
                             std::vector<double>& flux) const
{
    auto v = std::vector<double>(this->nnc_.numConnections(), 0.0);

    for (const auto& cat : this->nnc_.allCategories()) {
        const auto& rel    = this->nnc_.getRelations(cat);
        const auto  fluxID = vector + rel.fluxID;

        auto gridID = 0;
        for (const auto& G : this->grid_) {
            const auto& iset = rel.indexSet.getGridCollection(gridID);

            // Must increment grid ID irrespective of early break of
            // iteration.
            gridID += 1;

            if (iset.empty()) {
                // No NNCs for this category in this grid.  Skip.
                continue;
            }

            // Note: Method name is confusing, but does actually do what we
            // want here.
            const auto q = G.cellData(this->src_.get(), fluxID);

            if (q.empty()) {
                // No flux data for this category in this grid.  Skip.
                continue;
            }

            // Data fully available for (category,gridID).  Assign
            // approriate subset of NNC flux vector.
            for (const auto& ix : iset) {
                assert (ix.neighIdx < v.size());
                assert (ix.kwIdx    < q.size());

                v[ix.neighIdx] = q[ix.kwIdx];
            }
        }
    }

    flux.insert(flux.end(), v.begin(), v.end());
}

std::string
Opm::ECLGraph::Impl::
flowVector(const BlackoilPhases::PhaseIndex phase) const
{
    const auto vector = std::string("FLR"); // Flow-rate, reservoir

    if (phase == BlackoilPhases::PhaseIndex::Aqua) {
        return vector + "WAT";
    }

    if (phase == BlackoilPhases::PhaseIndex::Liquid) {
        return vector + "OIL";
    }

    if (phase == BlackoilPhases::PhaseIndex::Vapour) {
        return vector + "GAS";
    }

    {
        std::ostringstream os;

        os << "Invalid phase index '" << phase << '\'';

        throw std::invalid_argument(os.str());
    }
}

// ======================================================================

Opm::ECLGraph::ECLGraph(ImplPtr pImpl)
    : pImpl_(std::move(pImpl))
{
}

Opm::ECLGraph::ECLGraph(ECLGraph&& rhs)
    : pImpl_(std::move(rhs.pImpl_))
{}

Opm::ECLGraph::~ECLGraph()
{}

Opm::ECLGraph&
Opm::ECLGraph::operator=(ECLGraph&& rhs)
{
    this->pImpl_ = std::move(rhs.pImpl_);

    return *this;
}

Opm::ECLGraph
Opm::ECLGraph::load(const Path& grid, const Path& init)
{
    auto pImpl = ImplPtr{new Impl(grid, init)};

    return { std::move(pImpl) };
}

void
Opm::ECLGraph::assignFluxDataSource(const Path& src)
{
    this->pImpl_->assignDataSource(src);
}

std::size_t
Opm::ECLGraph::numCells() const
{
    return this->pImpl_->numCells();
}

std::size_t
Opm::ECLGraph::numConnections() const
{
    return this->pImpl_->numConnections();
}

std::vector<int> Opm::ECLGraph::neighbours() const
{
    return this->pImpl_->neighbours();
}

std::vector<double> Opm::ECLGraph::poreVolume() const
{
    return this->pImpl_->activePoreVolume();
}

std::vector<double>
Opm::ECLGraph::
flux(const BlackoilPhases::PhaseIndex phase,
     const int                        rptstep) const
{
    return this->pImpl_->flux(phase, rptstep);
}
