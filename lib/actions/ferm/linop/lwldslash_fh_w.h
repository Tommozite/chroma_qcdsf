// -*- C++ -*-
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_fh_h__
#define __lwldslash_fh_h__

#include "state.h"
#include "io/aniso_io.h"
#include "actions/ferm/linop/lwldslash_base_w.h"

namespace Chroma
{
    //! General Wilson-Dirac dslash
    /*!
   * \ingroup linop
   *
   * DSLASH
   *
   * This routine is specific to Wilson fermions!
   *
   * Description:
   *
   * This routine applies the operator D' to Psi, putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \
   *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
   *	       /    mu			  mu
   *	       ---
   *	       mu=0
   *
   *	             Nd-1
   *	             ---
   *	             \    +
   *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *	             /    mu			   mu
   *	             ---
   *	             mu=0
   *
   */

    template <typename T, typename P, typename Q>
    class QDPWilsonDslashFHT : public WilsonDslashBase<T, P, Q>
    {
    public:
        //! Empty constructor. Must use create later
        QDPWilsonDslashFHT();

        //! Full constructor
        QDPWilsonDslashFHT(Handle<FermState<T, P, Q>> state);

        //! Full constructor with anisotropy
        QDPWilsonDslashFHT(Handle<FermState<T, P, Q>> state,
                           const AnisoParam_t &aniso_);

        //! Full constructor with general coefficients
        QDPWilsonDslashFHT(Handle<FermState<T, P, Q>> state,
                           const multi1d<Real> &coeffs_);

        //! Creation routine
        void create(Handle<FermState<T, P, Q>> state);

        //! Creation routine with anisotropy
        void create(Handle<FermState<T, P, Q>> state,
                    const AnisoParam_t &aniso_);

        //! Full constructor with general coefficients
        void create(Handle<FermState<T, P, Q>> state,
                    const multi1d<Real> &coeffs_);

        //! No real need for cleanup here
        ~QDPWilsonDslashFHT() {}

        /**
     * Apply a dslash
     *
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
     *
     * \return The output of applying dslash on psi
     */
        void apply(T &chi, const T &psi, enum PlusMinus isign, int cb) const;

        //! Return the fermion BC object for this linear operator
        const FermBC<T, P, Q> &getFermBC() const { return *fbc; }

    protected:
        //! Get the anisotropy parameters
        const multi1d<Real> &getCoeffs() const { return coeffs; }

    private:
        multi1d<Real> coeffs; /*!< Nd array of coefficients of terms in the action */
        Handle<FermBC<T, P, Q>> fbc;
        Q u;
    };

    //! General Wilson-Dirac dslash
    /*! \ingroup linop
   * DSLASH
   *
   * This routine is specific to Wilson fermions!
   *
   * Description:
   *
   * This routine applies the operator D' to Psi, putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \
   *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
   *	       /    mu			  mu
   *	       ---
   *	       mu=0
   *
   *	             Nd-1
   *	             ---
   *	             \    +
   *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *	             /    mu			   mu
   *	             ---
   *	             mu=0
   *
   */

    //! Empty constructor
    template <typename T, typename P, typename Q>
    QDPWilsonDslashFHT<T, P, Q>::QDPWilsonDslashFHT() {}

    //! Full constructor
    template <typename T, typename P, typename Q>
    QDPWilsonDslashFHT<T, P, Q>::QDPWilsonDslashFHT(Handle<FermState<T, P, Q>> state)
    {
        create(state);
    }

    //! Full constructor with anisotropy
    template <typename T, typename P, typename Q>
    QDPWilsonDslashFHT<T, P, Q>::QDPWilsonDslashFHT(Handle<FermState<T, P, Q>> state,
                                                    const AnisoParam_t &aniso_)
    {
        create(state, aniso_);
    }

    //! Full constructor with general coefficients
    template <typename T, typename P, typename Q>
    QDPWilsonDslashFHT<T, P, Q>::QDPWilsonDslashFHT(Handle<FermState<T, P, Q>> state,
                                                    const multi1d<Real> &coeffs_)
    {
        create(state, coeffs_);
    }

    //! Creation routine
    template <typename T, typename P, typename Q>
    void QDPWilsonDslashFHT<T, P, Q>::create(Handle<FermState<T, P, Q>> state)
    {
        multi1d<Real> cf(Nd);
        cf = 1.0;
        create(state, cf);
    }

    //! Creation routine with anisotropy
    template <typename T, typename P, typename Q>
    void QDPWilsonDslashFHT<T, P, Q>::create(Handle<FermState<T, P, Q>> state,
                                             const AnisoParam_t &anisoParam)
    {
        START_CODE();

        create(state, makeFermCoeffs(anisoParam));

        END_CODE();
    }

    //! Full constructor with general coefficients
    template <typename T, typename P, typename Q>
    void QDPWilsonDslashFHT<T, P, Q>::create(Handle<FermState<T, P, Q>> state,
                                             const multi1d<Real> &coeffs_)
    {
        //QDPIO::cout << "Setting up QDP Wilson Dslash\n";

        // Save a copy of the aniso params original fields and with aniso folded in
        coeffs = coeffs_;

        // Save a copy of the fermbc
        fbc = state->getFermBC();

        // Sanity check
        if (fbc.operator->() == 0)
        {
            QDPIO::cerr << "QDPWilsonDslashFH: error: fbc is null" << std::endl;
            QDP_abort(1);
        }

        u.resize(Nd);

        // Fold in anisotropy
        for (int mu = 0; mu < u.size(); ++mu)
        {
            u[mu] = (state->getLinks())[mu];
        }

        // Rescale the u fields by the anisotropy
        for (int mu = 0; mu < u.size(); ++mu)
        {
            u[mu] *= coeffs[mu];
        }
    }

    //! General Wilson-Dirac dslash
    /*! \ingroup linop
   * Wilson dslash
   *
   * Arguments:
   *
   *  \param chi	      Result				                (Write)
   *  \param psi	      Pseudofermion field				(Read)
   *  \param isign      D'^dag or D' ( MINUS | PLUS ) resp.		(Read)
   *  \param cb	      Checkerboard of OUTPUT std::vector			(Read) 
   */
    template <typename T, typename P, typename Q>
    void
    QDPWilsonDslashFHT<T, P, Q>::apply(T &chi, const T &psi,
                                       enum PlusMinus isign, int cb) const
    {
        START_CODE();
#if (QDP_NC == 2) || (QDP_NC == 3)
        /*     F 
     *   a2  (x)  :=  U  (x) (1 - isign gamma  ) psi(x)
     *     mu          mu                    mu
     */
        /*     B           +
     *   a2  (x)  :=  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
     *     mu          mu                       mu
     */
        // Recontruct the bottom two spinor components from the top two
        /*                        F           B
     *   chi(x) :=  sum_mu  a2  (x)  +  a2  (x)
     *                        mu          mu
     */

        /* Why are these lines split? An array syntax would help, but the problem is deeper.
     * The expression templates require NO variable args (like int's) to a function
     * and all args must be known at compile time. Hence, the function names carry
     * (as functions usually do) the meaning (and implicit args) to a function.
     */
        switch (isign)
        {
        case PLUS:
            chi[rb[cb]] = u[0] * shift(spinProjectDir0MinusFull(psi), FORWARD, 0) + shift(adj(u[0]) * spinProjectDir0PlusFull(psi), BACKWARD, 0)
#if QDP_ND >= 2
                          + u[1] * shift(spinProjectDir1MinusFull(psi), FORWARD, 1) + shift(adj(u[1]) * spinProjectDir1PlusFull(psi), BACKWARD, 1)
#endif
#if QDP_ND >= 3
                          + u[2] * shift(spinProjectDir2MinusFull(psi), FORWARD, 2) + shift(adj(u[2]) * spinProjectDir2PlusFull(psi), BACKWARD, 2)
#endif
#if QDP_ND >= 4
                          + u[3] * shift(spinProjectDir3MinusFull(psi), FORWARD, 3) + shift(adj(u[3]) * spinProjectDir3PlusFull(psi), BACKWARD, 3)
#endif

                //             chi[rb[cb]] = spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(psi), FORWARD, 0)) + spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(psi), BACKWARD, 0))
                // #if QDP_ND >= 2
                //                           + spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(psi), FORWARD, 1)) + spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(psi), BACKWARD, 1))
                // #endif
                // #if QDP_ND >= 3
                //                           + spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(psi), FORWARD, 2)) + spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(psi), BACKWARD, 2))
                // #endif
                // #if QDP_ND >= 4
                //                           + spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(psi), FORWARD, 3)) + spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(psi), BACKWARD, 3))
                // #endif
                // #if QDP_ND >= 5
                // #error "Unsupported number of dimensions"
                // #endif
                ;
            break;

        case MINUS:
            chi[rb[cb]] = u[0] * shift(spinProjectDir0PlusFull(psi), FORWARD, 0) + shift(adj(u[0]) * spinProjectDir0MinusFull(psi), BACKWARD, 0)
#if QDP_ND >= 2
                          + u[1] * shift(spinProjectDir1PlusFull(psi), FORWARD, 1) + shift(adj(u[1]) * spinProjectDir1MinusFull(psi), BACKWARD, 1)
#endif
#if QDP_ND >= 3
                          + u[2] * shift(spinProjectDir2PlusFull(psi), FORWARD, 2) + shift(adj(u[2]) * spinProjectDir2MinusFull(psi), BACKWARD, 2)
#endif
#if QDP_ND >= 4
                          + u[3] * shift(spinProjectDir3PlusFull(psi), FORWARD, 3) + shift(adj(u[3]) * spinProjectDir3MinusFull(psi), BACKWARD, 3)
#endif
//             chi[rb[cb]] = spinReconstructDir0Plus(u[0] * shift(spinProjectDir0Plus(psi), FORWARD, 0)) + spinReconstructDir0Minus(shift(adj(u[0]) * spinProjectDir0Minus(psi), BACKWARD, 0))
// #if QDP_ND >= 2
//                           + spinReconstructDir1Plus(u[1] * shift(spinProjectDir1Plus(psi), FORWARD, 1)) + spinReconstructDir1Minus(shift(adj(u[1]) * spinProjectDir1Minus(psi), BACKWARD, 1))
// #endif
// #if QDP_ND >= 3
//                           + spinReconstructDir2Plus(u[2] * shift(spinProjectDir2Plus(psi), FORWARD, 2)) + spinReconstructDir2Minus(shift(adj(u[2]) * spinProjectDir2Minus(psi), BACKWARD, 2))
// #endif
// #if QDP_ND >= 4
//                           + spinReconstructDir3Plus(u[3] * shift(spinProjectDir3Plus(psi), FORWARD, 3)) + spinReconstructDir3Minus(shift(adj(u[3]) * spinProjectDir3Minus(psi), BACKWARD, 3))
// #endif
#if QDP_ND >= 5
#error "Unsupported number of dimensions"
#endif
                ;
            break;
        }

        QDPWilsonDslashFHT<T, P, Q>::getFermBC().modifyF(chi, QDP::rb[cb]);
#else
        QDPIO::cerr
            << "lwldslash_w_fh: not implemented for NC!=3\n";
        QDP_abort(13);
#endif
        END_CODE();
    }

    typedef QDPWilsonDslashFHT<LatticeFermion,
                               multi1d<LatticeColorMatrix>,
                               multi1d<LatticeColorMatrix>>
        QDPWilsonDslashFH;

    typedef QDPWilsonDslashFHT<LatticeFermionF,
                               multi1d<LatticeColorMatrixF>,
                               multi1d<LatticeColorMatrixF>>
        QDPWilsonDslashFHF;

    typedef QDPWilsonDslashFHT<LatticeFermionD,
                               multi1d<LatticeColorMatrixD>,
                               multi1d<LatticeColorMatrixD>>
        QDPWilsonDslashFHD;

} // End Namespace Chroma

#endif
