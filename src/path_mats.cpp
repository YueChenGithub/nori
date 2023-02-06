#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

    class PathMats : public Integrator {
    public:
        PathMats(const PropertyList& props) {}

        Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& _ray) const override {
            Color3f color = 0;//最终结果
            Color3f t = 1;//本次弹射对最终结果的贡献
            Ray3f rayRecursive = _ray;
            float probability;//本次继续的概率
            int depth = 1;//递归深度
            while (true) {
                Intersection its;
                if (!scene->rayIntersect(rayRecursive, its))
                    break;
                //光源贡献
                if (its.mesh->isEmitter()) {
                    EmitterQueryRecord lRecE(rayRecursive.o, its.p, its.shFrame.n);
                    color += t * its.mesh->getEmitter()->eval(lRecE);
                }
                //俄罗斯轮盘赌
                if (depth >= 3) {
                    probability = std::min(t.maxCoeff(), 0.99f);//选择路径贡献中最大项
                    if (sampler->next1D() > probability)
                        break;
                    t /= probability;
                }
                BSDFQueryRecord bRec(its.shFrame.toLocal(-rayRecursive.d));
                Color3f f = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                t *= f;//将贡献叠加到路径上
                //继续递归
                rayRecursive = Ray3f(its.p, its.toWorld(bRec.wo));
                depth++;
            }
            return color;
        }

        std::string toString() const {
            return "PathMats[]";
        }
    };

    NORI_REGISTER_CLASS(PathMats, "path_mats");
NORI_NAMESPACE_END