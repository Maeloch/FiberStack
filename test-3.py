import pybullet as p
import numpy as np
c = p.connect(p.DIRECT)
p.setGravity(0, 0, -9.81e-3)
p.setPhysicsEngineParameter(fixedTimeStep=2e-3, numSubSteps=4, numSolverIterations=50)
col = p.createCollisionShape(p.GEOM_CAPSULE, radius=0.01, height=1.0)
b = p.createMultiBody(baseMass=1.6e-4, baseCollisionShapeIndex=col,
                      basePosition=[0, 0, 1.0])
p.changeDynamics(b, -1, linearDamping=0.0023, angularDamping=0.0023,
                 restitution=0.0)
p.resetBaseVelocity(b, [0.05, 0, 0.0], [0, 0, 0])   # dérive lente, style rebond résiduel
for k in range(6000):
    p.stepSimulation()
    if k % 500 == 0:
        pos, _ = p.getBasePositionAndOrientation(b)
        lin, _ = p.getBaseVelocity(b)
        print(f"pas {k:5d}  z={pos[2]:.4f} mm  v={np.linalg.norm(lin):.4f} mm/ms")
