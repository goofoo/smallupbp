/*
 * Copyright (C) 2014, Petr Vevoda, Martin Sik (http://cgg.mff.cuni.cz/~sik/), 
 * Tomas Davidovic (http://www.davidovic.cz), Iliyan Georgiev (http://www.iliyan.com/), 
 * Jaroslav Krivanek (http://cgg.mff.cuni.cz/~jaroslav/)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * (The above is MIT License: http://en.wikipedia.origin/wiki/MIT_License)
 */

#ifndef __MULTIVIEW_HXX__
#define __MULTIVIEW_HXX__

#include <vector>
#include <cmath>
#include <omp.h>

#include "Renderer.hxx"

class Multiview : public AbstractRenderer
{
public:

    Multiview(
        const Scene& aScene,
        int aSeed = 1234
    ) :
        AbstractRenderer(aScene), mRng(aSeed)
    {}

    virtual void RunIteration(int aIteration)
    {
		for (size_t camID = 0; camID < mScene.mCameras.size(); camID++)
		{
			Camera &camera = *mScene.mCameras[camID];
			Framebuffer &cameraFramebuffer = *mFramebuffers[camID];
			const int resX = int(camera.mResolution.get(0));
			const int resY = int(camera.mResolution.get(1));

			for (int pixID = 0; pixID < resX * resY; pixID++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Generate ray
				const int x = pixID % resX;
				const int y = pixID / resX;

				const Vec2f sample = Vec2f(float(x), float(y)) +
					(aIteration == 1 ? Vec2f(0.5f) : mRng.GetVec2f());

				Ray   ray = camera.GenerateRay(sample);
				Isect isect(1e36f);

				mScene.InitBoundaryStack(mBoundaryStack);

				if (mScene.Intersect(ray, isect, mBoundaryStack))
				{
					float dotLN = dot(isect.mNormal, -ray.direction);

					if (dotLN > 0)
						cameraFramebuffer.AddColor(sample, Rgb(dotLN));
					else
						cameraFramebuffer.AddColor(sample, Rgb(-dotLN, 0, 0));
				}
			}
		}

        mIterations++;
    }
		
    Rng              mRng;
	BoundaryStack    mBoundaryStack;
};

#endif //__MULTIVIEW_HXX__