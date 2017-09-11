/*
 * MIT License: http://en.wikipedia.origin/wiki/MIT_License
 */

#ifndef __MULTIVIEW_HXX__
#define __MULTIVIEW_HXX__

#include <vector>
#include <cmath>
#include <omp.h>

#include "..\Beams\PhBeams.hxx"
#include "..\Misc\Timer.hxx"
#include "..\Path\VltPathVertex.hxx"
#include "..\Path\StaticArray.hxx"
#include "Renderer.hxx"



class Multiview : public AbstractRenderer
{
	// The sole point of this structure is to make carrying around the ray baggage easier.
	struct SubPathState
	{
		Pos   mOrigin;            // Path origin
		Dir   mDirection;         // Where to go next
		Rgb   mThroughput;        // Path throughput
		uint  mPathLength;        // Number of path segments, including this
		uint  mIsFiniteLight;     // Just generate by finite light
		uint  mSpecularPath : 1;  // All scattering events so far were specular

		Isect mIsect;
		PhotonBeamsArray mBeams;
		BoundaryStack mBoundaryStack; // Boundary stack
		VolumeSegments mVolumnSegments;
		LiteVolumeSegments mLiteVolumnSegments;
	};

	typedef std::vector<SubPathState> FullPath;

public:

	enum Version
	{
		kLightTracer = 0x0001,
		kPointPoint3D = 0x0010,
		kPointPoint = kPointPoint3D,
		kPointBeam2D = 0x0100,
		kPointBeam3D = 0x0200,
		kPointBeam = kPointBeam2D | kPointBeam3D,
		kBeamBeam1D = 0x1000,
		kBeamBeam2D = 0x2000,
		kBeamBeam3D = 0x4000,
		kBeamBeam = kBeamBeam1D | kBeamBeam2D | kBeamBeam3D,
		kPhotons = kPointPoint | kPointBeam | kBeamBeam,
		kAll = kLightTracer | kBeamBeam
	};

	Multiview(
		const Scene&            aScene,
		const int               aSeed,
		const Version           aVersion,

		const float			    aBB1DRadiusInitial,
		const float			    aBB1DRadiusAlpha,
		const RadiusCalculation aBB1DRadiusCalculation,
		const int  			    aBB1DRadiusKNN,
		const BeamType		    aBB1DBeamType,
		const float             aBB1DUsedLightSubPathCount,

		const float             aRefPathCountPerIter,
		const bool              aVerbose = true
	)
		: AbstractRenderer(aScene)
		, mRng(aSeed)
		, mVersion(aVersion)
		, mGlobalMediumID(aScene.mGlobalMediumID)
		, mPhotonBeams(aScene)
		, mVerbose(aVerbose)
		, mBB1DRadiusInitial(aBB1DRadiusInitial)
		, mBB1DRadiusAlpha(aBB1DRadiusAlpha)
		, mBB1DRadiusCalculation(aBB1DRadiusCalculation)
		, mBB1DRadiusKNN(aBB1DRadiusKNN)
		, mBB1DBeamType(aBB1DBeamType)
		, mBB1DUsedLightSubPathCount(aBB1DUsedLightSubPathCount)
		, mRefPathCountPerIter(aRefPathCountPerIter)
	{
		if (mBB1DRadiusInitial < 0)
			mBB1DRadiusInitial = -mBB1DRadiusInitial * mScene.mSceneSphere.mSceneRadius;
		UPBP_ASSERT(mBB1DRadiusInitial > 0);

		const Camera &tmpCamera = *mScene.mCameras[0];
		mLightSubPathCount = tmpCamera.mResolution.get(0) * tmpCamera.mResolution.get(1);

		if (mBB1DUsedLightSubPathCount < 0)
			mBB1DUsedLightSubPathCount = -mBB1DUsedLightSubPathCount * mLightSubPathCount;;

		// Reserve space
		mPhotonBeamsArray.reserve(mLightSubPathCount * 2);
		mFullPathes.reserve(mRefPathCountPerIter * 2);

		mPhotonBeams.mSeed = aSeed;
		for (size_t camId = 0; camId < mScene.mCameras.size(); camId++)
		{
			mAccumCameraFramebuffers.push_back(new Framebuffer());
			Vec2f resolution = mScene.mCameras[camId]->mResolution;
			size_t NPixels = size_t(resolution.x * resolution.y);
			mCameraPhotonMasks.push_back(std::vector<int>(NPixels, 0));
		}
	}

	~Multiview()
	{
		for (size_t camId = 0; camId < mScene.mCameras.size(); camId++)
		{
			delete mAccumCameraFramebuffers[camId];
		}
	}

	virtual void RunIteration(int aIteration)
	{
		Timer lightSubpathTimer;
		if (mVerbose)
			std::cout << " + tracing light sub-paths..." << std::endl;
		lightSubpathTimer.Start();

		// Radius reduction (1st iteration has aIteration == 0, thus offset)
		float radiusBB1D = mBB1DRadiusInitial * std::pow(1 + aIteration * mBB1DUsedLightSubPathCount / mRefPathCountPerIter, mBB1DRadiusAlpha - 1);
		radiusBB1D = std::max(radiusBB1D, 1e-7f); // Purely for numeric stability

		//////////////////////////////////////////////////////////////////////////
		// Generate light paths
		//////////////////////////////////////////////////////////////////////////
		int ABCSurfMutate = 0;
		int ABCMediaMutate = 0;
		int ABMutate = 0;
		int WillMutate = 0;
		if (mScene.GetLightCount() > 0 && mMaxPathLength > 1)
		{
			// If this renderer(thread) first be called
			if (mIterations == 0)
			{
				for (int pathIdx = 0; pathIdx < mLightSubPathCount; pathIdx++)
				{
					FullPath fullpath;
					mFullPathes.push_back(fullpath);

					// Generate light path origin and direction
					SubPathState lightState;
					GenerateLightSample(lightState);
					
					// In attenuating media the ray can never travel from infinity
					if (!lightState.mIsFiniteLight && mScene.GetGlobalMediumPtr()->HasAttenuation())
						continue;

					// We assume that the light is on surface
					bool originInMedium = false;

					//////////////////////////////////////////////////////////////////////////
					// Trace light path
					for (;; ++lightState.mPathLength)
					{
						// Prepare ray
						Ray ray(lightState.mOrigin, lightState.mDirection);
						Isect isect(1e36f);
						
						// Trace ray
						VolumeSegments subpathVolumeSegments;
						LiteVolumeSegments subpathLiteVolumeSegments;
						subpathVolumeSegments.clear();
						subpathLiteVolumeSegments.clear();

						bool intersected = mScene.Intersect(ray, originInMedium ? AbstractMedium::kOriginInMedium : 0, mRng, isect, lightState.mBoundaryStack, subpathVolumeSegments, subpathLiteVolumeSegments);

						lightState.mIsect = isect;
						CloneVolSegs(lightState.mVolumnSegments, subpathVolumeSegments);
						CloneLiteVolSegs(lightState.mLiteVolumnSegments, subpathLiteVolumeSegments);

						mFullPathes[pathIdx].push_back(lightState);

						// Store beam if required
						if ((mVersion & kBeamBeam) && pathIdx < mBB1DUsedLightSubPathCount)
						{
							AddBeams(ray, mFullPathes[pathIdx].back(), subpathVolumeSegments, subpathLiteVolumeSegments);
						}

						if (!intersected)
							break;

						UPBP_ASSERT(isect.IsValid());

						// Attenuate by intersected media (if any)
						if (!subpathVolumeSegments.empty())
						{
							// PDF
							float raySamplePdf = VolumeSegment::AccumulatePdf(subpathVolumeSegments);
							UPBP_ASSERT(raySamplePdf > 0);

							// Attenuation
							lightState.mThroughput *= VolumeSegment::AccumulateAttenuationWithoutPdf(subpathVolumeSegments) / raySamplePdf;
						}

						if (lightState.mThroughput.isBlackOrNegative())
							break;

						// Prepare scattering function at the hitpoint (BSDF/phase depending on whether the hitpoint is at surface or in media, the isect knows)
						BSDF bsdf(ray, isect, mScene, BSDF::kFromLight, mScene.RelativeIOR(isect, lightState.mBoundaryStack));

						if (!bsdf.IsValid()) // e.g. hitting surface too parallel with tangent plane
							break;

						// Compute hitpoint
						const Pos hitPoint = ray.origin + ray.direction * isect.mDist;

						// Terminate if the path would become too long after scattering
						if (lightState.mPathLength + 2 > mMaxPathLength)
							break;

						// Continue random walk
						if (!SampleScattering(bsdf, hitPoint, lightState, isect))
							break;
					} // Trace light path end
					  ////////////////////////////////////////////////////////////////////////////
				}
			}
			else
			{
				// mIterations > 0
				for (size_t fpIdx = 0; fpIdx < mFullPathes.size(); fpIdx++)
				{
					FullPath &pickedFullPath = mFullPathes[fpIdx];
					int mutateIdx = mRng.GetFloat() * float(pickedFullPath.size());
					bool pointAInMedium = false;
					/*
					*   A              A
					*    \     C        \     C
					*     \   /    =>    \   /
					*      \ /            \ /
					*       B              Z
					*/
					if (mutateIdx == 0)
					{
						pointAInMedium = false; // assume light is not in medium
					}
					else
					{
						pointAInMedium = !pickedFullPath[mutateIdx - 1].mIsect.IsOnSurface();
					}

					SubPathState &subPathAB = pickedFullPath[mutateIdx];
					SubPathState subPathAZ;
					CloneSubPathState(subPathAZ, subPathAB);
					subPathAZ.mBeams.clear();

					bool willMutate = TryGetMutatedSubpath(subPathAZ, subPathAB, pointAInMedium);
					if (willMutate)
					{
						WillMutate++;
						if (mutateIdx + 1 < pickedFullPath.size()) // BC path Exists
						{
							Ray rayAZ(subPathAZ.mOrigin, subPathAZ.mDirection);
							Pos hitPointZ = rayAZ.target(subPathAZ.mIsect.mDist);
							BSDF bsdfZ(rayAZ, subPathAZ.mIsect, mScene, BSDF::kFromLight, mScene.RelativeIOR(subPathAZ.mIsect, subPathAZ.mBoundaryStack));
							
							SubPathState &subPathBC = pickedFullPath[mutateIdx + 1];
							if (bsdfZ.IsValid() && subPathBC.mIsect.IsValid()) //C point Exists
							{
								Pos hitPointC = subPathBC.mOrigin + subPathBC.mIsect.mDist * subPathBC.mDirection;
								Dir dirZC = hitPointC - hitPointZ;
								Ray rayZC(hitPointZ, dirZC.getNormalized());

								Rgb throughputZ = subPathAZ.mThroughput;
								AttenuateByIntersectedMedia(throughputZ, subPathAZ.mVolumnSegments);

								if (bsdfZ.IsOnSurface() && !throughputZ.isBlackOrNegative())
								{	
									// mutate scattering
									float cosThetaOut, sfDirPdfW;
									Rgb sfFactor = bsdfZ.Evaluate(rayZC.direction, cosThetaOut, &sfDirPdfW);
									throughputZ *= sfFactor * (cosThetaOut / sfDirPdfW);
									
									// only support diffuse
									// only mutate when Z's throughput min component > B's
									const float BETTER_THROUGHPUT_PROB = 0.7f;
									bool isToMutationbyBetterThroughput = false;
									if (throughputZ.isValid() && throughputZ.min() > subPathBC.mThroughput.min())
									{
										isToMutationbyBetterThroughput = mRng.GetFloat() < BETTER_THROUGHPUT_PROB;
									}

									if (throughputZ.isValid() && !throughputZ.isBlackOrNegative() && isToMutationbyBetterThroughput)
									{
										SubPathState subPathZC;
										CloneSubPathState(subPathZC, subPathBC);
										subPathZC.mBeams.clear();
										subPathZC.mOrigin = hitPointZ;
										subPathZC.mThroughput = throughputZ;
										subPathZC.mDirection = rayZC.direction;
										Isect isectZC(dirZC.size());
										BoundaryStack boundaryStackZC;
										mScene.InitBoundaryStack(boundaryStackZC);
										VolumeSegments volumeSegmentsZC;
										LiteVolumeSegments liteVolumeSegmentsZC;

										bool intersected = mScene.Intersect(rayZC, 0, mRng, isectZC, boundaryStackZC, volumeSegmentsZC, liteVolumeSegmentsZC);
										if ((rayZC.target(isectZC.mDist) - hitPointC).size() < 1e-6)
										{
											subPathZC.mIsect = isectZC;
											subPathZC.mBoundaryStack = boundaryStackZC;
											CloneVolSegs(subPathZC.mVolumnSegments,volumeSegmentsZC);
											CloneLiteVolSegs(subPathZC.mLiteVolumnSegments, liteVolumeSegmentsZC);

											AddBeams(rayZC, subPathZC, subPathZC.mVolumnSegments, subPathZC.mLiteVolumnSegments);
											Rgb throughputC = subPathZC.mThroughput;
											AttenuateByIntersectedMedia(throughputC, subPathZC.mVolumnSegments);

											/* mutate success finally */
											ABCSurfMutate++;
											CloneSubPathState(subPathAB, subPathAZ);
											CloneSubPathState(subPathBC, subPathZC);
											UpdateMutateThroughput(pickedFullPath, mutateIdx + 1, throughputC);
										}
										else
										{
											/* Z don't reach C */
											// give up mutate
										}
									} // End !throughputZ.isBlackOrNegative()
								}
								else // in media
								{
									SubPathState subPathZC;
									CloneSubPathState(subPathZC, subPathBC);
									subPathZC.mBeams.clear();
									subPathZC.mOrigin = hitPointZ;
									subPathZC.mThroughput = throughputZ;
									subPathZC.mDirection = rayZC.direction;
									Isect isectZC(dirZC.size());
									BoundaryStack boundaryStackZC;
									mScene.InitBoundaryStack(boundaryStackZC);
									VolumeSegments volumeSegmentsZC;
									LiteVolumeSegments liteVolumeSegmentsZC;
									bool intersected = mScene.Intersect(rayZC, 0, mRng, isectZC, boundaryStackZC, volumeSegmentsZC, liteVolumeSegmentsZC);

									if ((rayZC.target(isectZC.mDist) - hitPointC).size() < 1e-6)
									{
										subPathZC.mIsect = isectZC;
										subPathZC.mBoundaryStack = boundaryStackZC;
										CloneVolSegs(subPathZC.mVolumnSegments,volumeSegmentsZC);
										CloneLiteVolSegs(subPathZC.mLiteVolumnSegments, liteVolumeSegmentsZC);

										AddBeams(rayZC, subPathZC, subPathZC.mVolumnSegments, subPathZC.mLiteVolumnSegments);

										/* if short beam, move B to Z and done */

										//Rgb throughputC = subPathZC.mThroughput;
										//if (!subPathZC.mVolumnSegments.empty())
										//{
										//	// PDF
										//	float raySamplePdf = VolumeSegment::AccumulatePdf(subPathZC.mVolumnSegments);
										//	// Attenuation
										//	throughputC *= VolumeSegment::AccumulateAttenuationWithoutPdf(subPathZC.mVolumnSegments) / raySamplePdf;
										//}

										/* mutate success finally */
										ABCMediaMutate++;
										CloneSubPathState(subPathAB, subPathAZ);
										CloneSubPathState(subPathBC, subPathZC);
									}
									else
									{
										/* Z don't reach C */
										// give up mutate
									}
								}
							}
						}
						else // mutSubPath has no next subpath
						{
							ABMutate++;
							CloneSubPathState(subPathAB, subPathAZ);
						}
					} // End willMutate
				}
			}
		}

		lightSubpathTimer.Stop();
		if (mVerbose)
			std::cout << "    - light sub-path tracing done in " << lightSubpathTimer.GetLastElapsedTime() << " sec. " << std::endl;
		std::cout << "WillMutate: " << WillMutate << std::endl;
		std::cout << "ABCSurfMutate: " << ABCSurfMutate << std::endl;
		std::cout << "ABCMediaMutate: " << ABCMediaMutate << std::endl;
		std::cout << "ABMutate: " << ABMutate << std::endl;

		if (mMaxPathLength > 1 && (mVersion & kBeamBeam))
		{
			CollectPhotonBeamsToArray();
			std::cout << "PhotonBeamsArraySize: " << mPhotonBeamsArray.size() << std::endl;
			mPhotonBeams.build(mPhotonBeamsArray, mBB1DRadiusCalculation, mBB1DRadiusInitial, mBB1DRadiusKNN, mVerbose);
		}

		//////////////////////////////////////////////////////////////////////////
		// Each camera tracing
		//////////////////////////////////////////////////////////////////////////

		for (size_t camId = 0; camId < mScene.mCameras.size(); camId++)
		{
			Timer cameraTimer;
			if (mVerbose)
				std::cout << " + camera[" << camId << "] tracing light sub-paths..." << std::endl;

			cameraTimer.Start();

			Camera &camera = *mScene.mCameras[camId];
			Framebuffer &cameraFramebuffer = *mFramebuffers[camId];

			// While we have the same number of pixels (camera paths)
			// and light paths, we do keep them separate for clarity reasons
			const int resX = int(camera.mResolution.get(0));
			const int resY = int(camera.mResolution.get(1));

			const int pathCount = resX * resY;

			//////////////////////////////////////////////////////////////////////////
			// Generate light paths
			//////////////////////////////////////////////////////////////////////////
			
			if (mScene.GetLightCount() > 0 && mMaxPathLength > 1)
			{
				for (int pathIdx = 0; pathIdx < mLightSubPathCount; pathIdx++)
				{
					// Generate light path origin and direction
					SubPathState lightState;
					GenerateLightSample(lightState);

					// In attenuating media the ray can never travel from infinity
					if (!lightState.mIsFiniteLight && mScene.GetGlobalMediumPtr()->HasAttenuation())
						continue;

					// We assume that the light is on surface
					bool originInMedium = false;

					//////////////////////////////////////////////////////////////////////////
					// Trace light path
					for (;; ++lightState.mPathLength)
					{
						// Prepare ray
						Ray ray(lightState.mOrigin, lightState.mDirection);
						Isect isect(1e36f);

						// Trace ray
						VolumeSegments subpathVolumeSegments;
						LiteVolumeSegments subpathLiteVolumeSegments;
						subpathVolumeSegments.clear();
						subpathLiteVolumeSegments.clear();

						bool intersected = mScene.Intersect(ray, originInMedium ? AbstractMedium::kOriginInMedium : 0, mRng, isect, lightState.mBoundaryStack, subpathVolumeSegments, subpathLiteVolumeSegments);

						if (!intersected)
							break;

						UPBP_ASSERT(isect.IsValid());

						// Same throughput before update by volume propagation
						const Rgb prevThroughput = lightState.mThroughput;

						// Attenuate by intersected media (if any)
						float raySamplePdf(1.0f);
						if (!subpathVolumeSegments.empty())
						{
							// PDF
							raySamplePdf = VolumeSegment::AccumulatePdf(subpathVolumeSegments);
							UPBP_ASSERT(raySamplePdf > 0);

							// Attenuation
							lightState.mThroughput *= VolumeSegment::AccumulateAttenuationWithoutPdf(subpathVolumeSegments) / raySamplePdf;
						}

						if (lightState.mThroughput.isBlackOrNegative())
							break;

						// Prepare scattering function at the hitpoint (BSDF/phase depending on whether the hitpoint is at surface or in media, the isect knows)
						BSDF bsdf(ray, isect, mScene, BSDF::kFromLight, mScene.RelativeIOR(isect, lightState.mBoundaryStack));

						if (!bsdf.IsValid()) // e.g. hitting surface too parallel with tangent plane
							break;

						// Compute hitpoint
						const Pos hitPoint = ray.origin + ray.direction * isect.mDist;

						// Connect to camera, unless scattering function is purely specular
						if ((mVersion & kLightTracer) && (!bsdf.IsDelta()))
						{
							if (lightState.mPathLength + 1 >= mMinPathLength)
								ConnectToCamera(camId, lightState, hitPoint, bsdf, isect);
						}

						// Terminate if the path would become too long after scattering
						if (lightState.mPathLength + 2 > mMaxPathLength)
							break;

						// Continue random walk
						if (!SampleScattering(bsdf, hitPoint, lightState, isect))
							break;
					}
				} // for each light sample from source
			}

			cameraTimer.Stop();
			if (mVerbose)
				std::cout << "    - camera[" << camId << "] light sub-path tracing done in " << cameraTimer.GetLastElapsedTime() << " sec. " << std::endl;

			//////////////////////////////////////////////////////////////////////////
			// Generate primary rays & accumulate "photon" contributions.
			//////////////////////////////////////////////////////////////////////////

			if ((mVersion & kBeamBeam) && !mPhotonBeamsArray.empty())
			{
				CalculatePhotonContributions(camId);
			}
		} // for each cameras

		if (mMaxPathLength > 1 && (mVersion & kBeamBeam))
		{
			mPhotonBeams.destroy();
		}

		mIterations++;
	}

private:

	//////////////////////////////////////////////////////////////////////////
	// Photon mapping methods
	//////////////////////////////////////////////////////////////////////////

	void CalculatePhotonContributions(const int aCameraID)
	{
		Timer photonTimer;
		Camera &camera = *mScene.mCameras[aCameraID];
		Framebuffer &cameraFramebuffer = *mFramebuffers[aCameraID];
		const int resX = int(camera.mResolution.get(0));
		const int resY = int(camera.mResolution.get(1));
		UPBP_ASSERT(mVersion & kPhotons);

		if (mVerbose)
			std::cout << " + camera[" << aCameraID << "] tracing primary rays..." << std::endl;
		photonTimer.Start();

// #pragma omp parallel for //useless
		for (int pixID = 0; pixID < resX * resY; pixID++)
		{
			Rng pixRng(pixID + aCameraID * resX * resY);
			const int x = pixID % resX;
			const int y = pixID / resX;

			// Generate pixel sample
			const Vec2f sample = Vec2f(float(x), float(y)) + pixRng.GetVec2f();

			// Create ray through the sample
			Ray ray = camera.GenerateRay(sample);
			Isect isect(1e36f);

			// Init the boundary stack with the global medium and add enclosing material and medium if present
			BoundaryStack pixBoundaryStack;
			mScene.InitBoundaryStack(pixBoundaryStack);
			if (camera.mMatID != -1 && camera.mMedID != -1)
				mScene.AddToBoundaryStack(camera.mMatID, camera.mMedID, pixBoundaryStack);

			bool  originInMedium = false;
			float rayLength;
			Rgb   surfaceRadiance(0.f), volumeRadiance(0.f);

			VolumeSegments pixVolumeSegments;
			pixVolumeSegments.clear();

			// Cast primary ray (do not scatter in media - pass all the way until we hit some object)
			bool hitSomething = mScene.Intersect(ray, isect, pixBoundaryStack, Scene::kSampleVolumeScattering, 0, &pixRng, &pixVolumeSegments, NULL, NULL);

			const AbstractMedium* medium = (isect.mMedID >= 0) ? mScene.GetMediumPtr(isect.mMedID) : nullptr;

			if (hitSomething)
			{
				rayLength = isect.mDist;
				if (isect.mLightID >= 0)
				{
					// Directly hit some light
					const AbstractLight *light = mScene.GetLightPtr(isect.mLightID);
					UPBP_ASSERT(light);

					// Compute directly emitted radiance
					float directIllumPdfA;
					surfaceRadiance = light->GetRadiance(mScene.mSceneSphere,
						ray.direction, ray.target(isect.mDist), &directIllumPdfA);

					// Attenuate surface radiance by media extinction
					if (medium)
					{
						surfaceRadiance *= medium->EvalAttenuation(ray, 0, isect.mDist);
					}

					UPBP_ASSERT(!surfaceRadiance.isNanInfNeg());
				}
			}
			else
			{
				// Didn't hit anything - we ignore background light for now.
				rayLength = INFINITY;
			}

			if ((mMaxPathLength > 1) && !pixVolumeSegments.empty())
			{
				for (size_t i = 0; i < pixVolumeSegments.size(); i++)
				{
					// -1 = no medium, 0 = clear medium
					if (pixVolumeSegments[i].mMediumID > 0) 
					{
						mCameraPhotonMasks[aCameraID][pixID] += 1;
					}
				}

				// UPBP_ASSERT(medium);

				if (mVersion & kBeamBeam)
				{
					GridStats gridStats;
					volumeRadiance = mPhotonBeams.evalBeamBeamEstimate(mBB1DBeamType, ray, pixVolumeSegments, BB1D, 0, NULL, &gridStats) / mBB1DUsedLightSubPathCount;

					mBeamDensity.Accumulate(pixID, gridStats);
				}

				UPBP_ASSERT(!volumeRadiance.isNanInfNeg());
			}

			cameraFramebuffer.AddColor(sample, surfaceRadiance + volumeRadiance);
		}

		photonTimer.Stop();
		if (mVerbose)
			std::cout << std::setprecision(3) << "   - camera[" << aCameraID << "] primary ray tracing done in " << photonTimer.GetLastElapsedTime() << " sec. " << std::endl;

		mCameraTracingTime += photonTimer.GetLastElapsedTime();
	}

	//////////////////////////////////////////////////////////////////////////
	// Light tracing methods
	//////////////////////////////////////////////////////////////////////////

	// Samples light emission
	void GenerateLightSample(SubPathState &oLightState)
	{
		// Sample lights uniformly

		const int   lightCount = mScene.GetLightCount();
		const float lightPickProb = 1.f / lightCount;

		const int   lightID = int(mRng.GetFloat() * lightCount);
		const AbstractLight *light = mScene.GetLightPtr(lightID);
		UPBP_ASSERT(light);

		// Generate light path origin and direction

		const Vec2f rndDirSamples = mRng.GetVec2f();
		const Vec2f rndPosSamples = mRng.GetVec2f();
		float emissionPdfW, directPdfW, cosLight;

		oLightState.mThroughput = light->Emit(mScene.mSceneSphere, rndDirSamples, rndPosSamples,
			oLightState.mOrigin, oLightState.mDirection,
			emissionPdfW, &directPdfW, &cosLight);

		// Complete light path state initialization

		emissionPdfW *= lightPickProb;
		directPdfW *= lightPickProb;  // (not used by pure light tracing)

		oLightState.mThroughput /= emissionPdfW;
		oLightState.mPathLength = 1;
		oLightState.mIsFiniteLight = light->IsFinite() ? 1 : 0;

		// Init the boundary stack with the global medium and add enclosing material and medium if present
		mScene.InitBoundaryStack(oLightState.mBoundaryStack);
		if (light->mMatID != -1 && light->mMedID != -1)
		{
			mScene.AddToBoundaryStack(light->mMatID, light->mMedID, oLightState.mBoundaryStack);
		}
	}

	// Computes contribution of light sample to camera by splatting is onto the
	// framebuffer. Multiplies by throughput (obviously, as nothing is returned).
	void ConnectToCamera(
		const int                        aCameraID,
		const SubPathState               &aLightState,
		const Pos                        &aHitpoint,
		const BSDF                       &aLightBSDF,
		const Isect                      &isect)
	{
		// Get camera and direction to it
		const Camera &camera = *mScene.mCameras[aCameraID];
		Framebuffer &cameraFramebuffer = *mFramebuffers[aCameraID];

		Dir directionToCamera = camera.mOrigin - aHitpoint;

		// Check point is in front of camera
		if (dot(camera.mDirection, -directionToCamera) <= 0.f)
			return;

		// Check it projects to the screen (and where)
		const Vec2f imagePos = camera.WorldToRaster(aHitpoint);
		if (!camera.CheckRaster(imagePos))
			return;

		// Compute distance and normalize direction to camera
		const float distEye2 = directionToCamera.square();
		const float distance = std::sqrt(distEye2);
		// UPBP_ASSERT((directionToCamera / distance) == directionToCamera.getNormalized());
		// directionToCamera /= distance;
		directionToCamera = directionToCamera.getNormalized();
		

		// Get the scattering function factor
		float cosToCamera, sfDirPdfW, sfRevPdfW;
		Rgb sfFactor = aLightBSDF.Evaluate(directionToCamera, cosToCamera, &sfDirPdfW, &sfRevPdfW);

		if (sfFactor.isBlackOrNegative())
			return;

		sfRevPdfW *= aLightBSDF.ContinuationProb();

		// Compute PDF conversion factor from image plane area to surface area
		const float cosAtCamera = dot(camera.mDirection, -directionToCamera);
		const float imagePointToCameraDist = camera.mImagePlaneDist / cosAtCamera;
		const float imageToSolidAngleFactor = Utils::sqr(imagePointToCameraDist) / cosAtCamera;
		const float imageToSurfaceFactor = imageToSolidAngleFactor * std::abs(cosToCamera) / Utils::sqr(distance);

		// We put the virtual image plane at such a distance from the camera origin
		// that the pixel area is one and thus the image plane sampling PDF is 1.
		// The area PDF of aHitpoint as sampled from the camera is then equal to
		// the conversion factor from image plane area density to surface area density
		const float cameraPdfA = imageToSurfaceFactor;

		const float surfaceToImageFactor = 1.f / imageToSurfaceFactor;

		VolumeSegments surfVolumeSegments;

		surfVolumeSegments.clear();
		if (!mScene.Occluded(aHitpoint, directionToCamera, distance, aLightState.mBoundaryStack, isect.IsOnSurface() ? 0 : AbstractMedium::kOriginInMedium, surfVolumeSegments))
		{
			// Get attenuation from intersected media (if any)
			Rgb mediaAttenuation = surfVolumeSegments.empty() ? Rgb(1.0) : VolumeSegment::AccumulateAttenuationWithoutPdf(surfVolumeSegments);

			// We divide the contribution by surfaceToImageFactor to convert the (already
			// divided) PDF from surface area to image plane area, w.r.t. which the
			// pixel integral is actually defined. We also divide by the number of samples
			// this technique makes, which is equal to the number of light sub-paths
			Rgb contrib = aLightState.mThroughput * sfFactor * mediaAttenuation / (mLightSubPathCount * surfaceToImageFactor);

			if (contrib.isBlackOrNegative())
				return;

			cameraFramebuffer.AddColor(imagePos, contrib);
		}
	}

	// Samples a scattering direction camera/light sample according to scattering function.
	// Returns false for termination
	bool SampleScattering(
		const BSDF      &aBSDF,
		const Pos       &aHitPoint,
		SubPathState    &aoState,
		const Isect     &isect)
	{
		// Sample scattering function		

		Dir   rndTriplet = mRng.GetVec3f(); // x,y for direction, z for component. No rescaling happens
		float sfDirPdfW, cosThetaOut;
		uint  sampledEvent;
		// BSDF.hxx#400 for media phase function
		// if isotrophic, PhaseFunction::Sample return RGB(1/4£k)
		Rgb sfFactor = aBSDF.Sample(rndTriplet, aoState.mDirection,
			sfDirPdfW, cosThetaOut, &sampledEvent);

		if (sfFactor.isBlackOrNegative())
			return false;

		bool specular = (sampledEvent & BSDF::kSpecular) != 0;

		// If we sampled specular event, then the reverse probability
		// cannot be evaluated, but we know it is exactly the same as
		// forward probability, so just set it. If non-specular event happened,
		// we evaluate the PDF
		float sfRevPdfW = sfDirPdfW;
		if (!specular)
			sfRevPdfW = aBSDF.Pdf(aoState.mDirection, BSDF::kReverse);

		// Russian roulette
		const float contProb = aBSDF.ContinuationProb();
		if (contProb == 0 || (contProb < 1.0f && mRng.GetFloat() > contProb))
			return false;

		sfDirPdfW *= contProb;
		sfRevPdfW *= contProb;

		// Update throughput
		aoState.mThroughput *= sfFactor * (cosThetaOut / sfDirPdfW);

		// Switch medium on refraction
		if ((sampledEvent & BSDF::kRefract) != 0)
			mScene.UpdateBoundaryStackOnRefract(isect, aoState.mBoundaryStack);

		// Update the rest of the light path state

		if (specular)
		{
			UPBP_ASSERT(sfDirPdfW == sfRevPdfW);
			aoState.mSpecularPath &= 1;
		}
		else
		{
			aoState.mSpecularPath &= 0;
		}

		aoState.mOrigin = aHitPoint;

		return true;
	}

	// Adds beams to beams array
	void AddBeams(
		const Ray		&aRay,
		SubPathState    &aLightState,
		const VolumeSegments &aVolumeSegments,
		const LiteVolumeSegments &aLiteVolumeSegments
	)
	{
		Rgb throughput = aLightState.mThroughput;
		if (mBB1DBeamType == SHORT_BEAM)
		{
			if (aVolumeSegments.size() > 0)
			{
				const VolumeSegment& first = aVolumeSegments[0];
				const VolumeSegment& cfirst = *(aVolumeSegments.cbegin());
				printf("");
			}
			
			for (VolumeSegments::const_iterator it = aVolumeSegments.cbegin(); it != aVolumeSegments.cend(); ++it)
			{
				UPBP_ASSERT(it->mMediumID >= 0);
				PhotonBeam beam;
				beam.mMedium = mScene.mMedia[it->mMediumID];
				if (beam.mMedium->HasScattering())
				{
					beam.mRay = Ray(aRay.origin + aRay.direction * it->mDistMin, aRay.direction);
					beam.mLength = it->mDistMax - it->mDistMin;
					beam.mFlags = SHORT_BEAM;
					beam.mThroughputAtOrigin = throughput;
					aLightState.mBeams.push_back(beam);
					// mPhotonBeamsArray.push_back(beam);
				}
				throughput *= it->mAttenuation / it->mRaySamplePdf;
				// short beam throughput *= (1, 1, 1)
			}
		}
		else // LONG_BEAM
		{
			UPBP_ASSERT(mBB1DBeamType == LONG_BEAM);
			for (LiteVolumeSegments::const_iterator it = aLiteVolumeSegments.cbegin(); it != aLiteVolumeSegments.cend(); ++it)
			{
				UPBP_ASSERT(it->mMediumID >= 0);
				PhotonBeam beam;
				beam.mMedium = mScene.mMedia[it->mMediumID];
				if (beam.mMedium->HasScattering())
				{
					beam.mRay = Ray(aRay.origin + aRay.direction * it->mDistMin, aRay.direction);
					beam.mLength = it->mDistMax - it->mDistMin;
					beam.mFlags = LONG_BEAM;
					beam.mThroughputAtOrigin = throughput;
					mPhotonBeamsArray.push_back(beam);
				}
				if (beam.mMedium->IsHomogeneous())
				{
					const HomogeneousMedium * medium = ((const HomogeneousMedium *)beam.mMedium);
					throughput *= medium->EvalAttenuation(it->mDistMax - it->mDistMin);
				}
				else
				{
					throughput *= beam.mMedium->EvalAttenuation(aRay, it->mDistMin, it->mDistMax);
				}
			}

			if (!mPhotonBeamsArray.empty() && mPhotonBeamsArray.back().mLength > mPhotonBeamsArray.back().mMedium->MaxBeamLength())
				// max(mPhotonBeamsArray.back().mLength, mPhotonBeamsArray.back().mMedium->MaxBeamLength())
				mPhotonBeamsArray.back().mLength = mPhotonBeamsArray.back().mMedium->MaxBeamLength();
		}
	}

	void CollectPhotonBeamsToArray()
	{
		mPhotonBeamsArray.clear();
		for (auto fp = mFullPathes.cbegin(); fp != mFullPathes.cend(); fp++)
		{
			const FullPath &fullpath = *fp;
			for (auto subpath = fullpath.cbegin(); subpath != fullpath.cend(); subpath++)
			{
				const PhotonBeamsArray &beams = subpath->mBeams;
				for (auto beam = beams.cbegin(); beam != beams.cend(); beam++)
				{
					mPhotonBeamsArray.push_back(*beam);
				}
			}
		}
	}

	bool IsSimilarSubPath(const SubPathState& a, const SubPathState& b)
	{
		// intersect on the same mat and the same med
		if (a.mIsect.mMatID != b.mIsect.mMatID || a.mIsect.mMedID != b.mIsect.mMedID)
		{
			return false;
		}			

		return true;
	}

	bool TryGetMutatedSubpath(SubPathState& mutatedSubpath, const SubPathState& aRefSubpath, bool aOriginInMedium) {
		
		int originInMedium = aOriginInMedium? AbstractMedium::kOriginInMedium : 0;

		if (aRefSubpath.mIsect.IsValid())
		{
			Pos hitPointB = aRefSubpath.mOrigin + aRefSubpath.mIsect.mDist * aRefSubpath.mDirection;
			Dir offsetBtoZ = (mRng.GetVec3f() - Dir(0.5)) * 0.1f; //¡Ó0.05
			mutatedSubpath.mDirection = (hitPointB + offsetBtoZ - aRefSubpath.mOrigin).getNormalized();

			Ray rayAZ(aRefSubpath.mOrigin, mutatedSubpath.mDirection);
			Isect isectAZ(1e36f);
			VolumeSegments volumeSegmentsAZ;
			LiteVolumeSegments liteVolumeSegmentsAZ;
			BoundaryStack boundaryStackAZ;
			mScene.InitBoundaryStack(boundaryStackAZ);

			bool intersected = mScene.Intersect(rayAZ, originInMedium, mRng, isectAZ, boundaryStackAZ, volumeSegmentsAZ, liteVolumeSegmentsAZ);

			mutatedSubpath.mIsect = isectAZ;
			mutatedSubpath.mBoundaryStack = boundaryStackAZ;
			CloneVolSegs(mutatedSubpath.mVolumnSegments, volumeSegmentsAZ);
			CloneLiteVolSegs(mutatedSubpath.mLiteVolumnSegments, liteVolumeSegmentsAZ);

			bool isSimilar = IsSimilarSubPath(mutatedSubpath, aRefSubpath);

			if (isSimilar)
			{
				if (mVersion & kBeamBeam)
				{
					AddBeams(rayAZ, mutatedSubpath, volumeSegmentsAZ, liteVolumeSegmentsAZ);
				}
			}

			return isSimilar;
		}
		else // !aRefSubpath.mIsect.IsValid()
		{
			return false;
		}
	}

	void ScatteringMutation(/*mutatedSubpath, relmutSubpath*/)
	{
	}

	void UpdateMutateThroughput(FullPath &aoPickedFullPath, int aIndexOfZC, const Rgb &aThroughputAfterZC)
	{
		if (aThroughputAfterZC.isBlackOrNegative())
		{
			// remove after
			if (aoPickedFullPath.begin() + aIndexOfZC != aoPickedFullPath.end())
			{
				aoPickedFullPath.erase(aoPickedFullPath.begin() + aIndexOfZC + 1, aoPickedFullPath.end());
			}
			return;
		}

		Rgb throughput(aThroughputAfterZC);

		for (size_t i = aIndexOfZC + 1; i < aoPickedFullPath.size(); i++)
		{
			SubPathState &sp = aoPickedFullPath[i];
			sp.mThroughput = throughput;
			sp.mBeams.clear();
			Ray ray(sp.mOrigin, sp.mDirection);
			AddBeams(ray, sp, sp.mVolumnSegments, sp.mLiteVolumnSegments);
			
			if (!sp.mVolumnSegments.empty())
			{
				// PDF
				float raySamplePdf = VolumeSegment::AccumulatePdf(sp.mVolumnSegments);
				// Attenuation
				throughput *= VolumeSegment::AccumulateAttenuationWithoutPdf(sp.mVolumnSegments) / raySamplePdf;
			}

			if (throughput.isBlackOrNegative())
			{
				// remove after
				aoPickedFullPath.erase(aoPickedFullPath.begin() + i + 1, aoPickedFullPath.end());
				break;
			}
		}
	}

	void AttenuateByIntersectedMedia(Rgb &aoThroughput, const VolumeSegments& aVolumeSegments)
	{
		if (!aVolumeSegments.empty())
		{
			// PDF
			float raySamplePdf = VolumeSegment::AccumulatePdf(aVolumeSegments);
			// Attenuation
			aoThroughput *= VolumeSegment::AccumulateAttenuationWithoutPdf(aVolumeSegments) / raySamplePdf;
		}
	}

	void CloneVolSegs(VolumeSegments &target,  const VolumeSegments &source)
	{
		target.clear();
		for (auto it = source.cbegin(); it != source.cend(); ++it)
		{
			target.push_back(*it);
		}
	}

	void CloneLiteVolSegs(LiteVolumeSegments &target, const LiteVolumeSegments &source)
	{
		target.clear();
		for (auto it = source.cbegin(); it != source.cend(); ++it)
		{
			target.push_back(*it);
		}
	}

	void CloneSubPathState(SubPathState &target, const SubPathState &source)
	{
		target = source;
		CloneVolSegs(target.mVolumnSegments, source.mVolumnSegments);
		CloneLiteVolSegs(target.mLiteVolumnSegments, source.mLiteVolumnSegments);
	}
private:
	// float mScreenPixelCount;    // Number of pixels // no use
	float mLightSubPathCount;   // Number of light sub-paths
	float mRefPathCountPerIter; // Reference number of paths per iteration

	std::vector<std::vector<int>> mCameraPhotonMasks;
	std::vector<Framebuffer *> mAccumCameraFramebuffers;
	std::vector<FullPath> mFullPathes;

	PhotonBeamsArray mPhotonBeamsArray;	    // Stores photon beams
	
	// VolumeSegments mVolumeSegments;         // Path segments intersecting media (up to scattering point)
	// LiteVolumeSegments mLiteVolumeSegments; // Lite path segments intersecting media (up to intersection with solid surface)

	PhotonBeamsEvaluator mPhotonBeams;

	int mGlobalMediumID;

	Rng               mRng;
	BoundaryStack     mBoundaryStack;
	Version           mVersion;
	bool              mVerbose;
	Timer             mTimer;

	// BB1D
	float                mBB1DRadiusInitial;         // Initial merging radius
	float                mBB1DRadiusAlpha;           // Radius reduction rate parameter
	RadiusCalculation    mBB1DRadiusCalculation;     // Type of photon radius calculation
	int	                 mBB1DRadiusKNN;             // Value x means that x-th closest beam vertex will be used for calculation of cone radius at the current beam vertex
	BeamType             mBB1DBeamType;              // Short/long beam
	float                mBB1DUsedLightSubPathCount; // First mBB1DUsedLightSubPathCount out of mLightSubPathCount light paths will generate photon beams
};

#endif //__MULTIVIEW_HXX__