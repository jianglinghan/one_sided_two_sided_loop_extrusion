import numpy as np

class LEFTranslocator:
    def __init__(
        self,
        numLEF,
        deathProb,
        stalledDeathProb,
        birthArray,
        pauseProb,
        stallProbLeft,
        stallProbRight,
        *args,
        mode="asymmetric",
        switchProb=0.0, # new params
        speedMultiplier=1.0,
        stability_factor=0.1,  
        dwell_time=10,
    ):
        self.numSite = len(birthArray)
        self.numLEF = numLEF
        self.mode = mode
        self.switchProb = switchProb
        self.speedMultiplier = speedMultiplier
        self.stability_factor = stability_factor 
        self.dwell_time = dwell_time

        self.stallProbLeft = stallProbLeft
        self.stallProbRight = stallProbRight
        self.deathProb = deathProb
        self.pause = pauseProb

        birthArray[0] = 0
        birthArray[-1] = 0
        self.birthProb = np.cumsum(birthArray, dtype=np.double)
        self.birthProb /= self.birthProb[-1]

        self.LEFs = np.zeros((self.numLEF, 2), int)
        self.stalled = np.zeros((self.numLEF, 2), int)
        self.directions = np.zeros(self.numLEF, int)
        self.currentSpeeds = np.ones(self.numLEF) 
        self.dwell = np.zeros(self.numLEF, int)

        self.occupied = np.zeros(self.numSite, int)
        self.stalledDeathProb = stalledDeathProb

        self.occupied[0] = 1
        self.occupied[-1] = 1

        for i in range(self.numLEF):
            self.LEF_birth(i)

    def LEF_birth(self, i):
        while True:
            pos = np.searchsorted(self.birthProb, np.random.random())
            if pos >= self.numSite - 1 or pos <= 0:
                continue
            if self.occupied[pos] == 1:
                continue
            self.LEFs[i] = pos
            self.occupied[pos] = 1
            self.dwell[i] = 0
            self.directions[i] = np.random.randint(0, 2)
            self.currentSpeeds[i] = 1.0 
            if (pos < (self.numSite - 3)) and (self.occupied[pos + 1] == 0):
                if np.random.random() > 0.5:
                    self.LEFs[i, 1] = pos + 1
                    self.occupied[pos + 1] = 1
            return

    def LEF_death(self):
        for i in range(self.numLEF):
            is_stable = (self.stalled[i, 0] == 1 and self.stalled[i, 1] == 1)
            
            if is_stable:
                deathProb = 0.0001 # or 0
            else:
                dp1 = self.deathProb[self.LEFs[i, 0]] if self.stalled[i, 0] == 0 else self.stalledDeathProb[self.LEFs[i, 0]]
                dp2 = self.deathProb[self.LEFs[i, 1]] if self.stalled[i, 1] == 0 else self.stalledDeathProb[self.LEFs[i, 1]]
                deathProb = max(dp1, dp2)
            if np.random.random() < deathProb:
                self.occupied[self.LEFs[i, 0]] = 0
                self.occupied[self.LEFs[i, 1]] = 0
                self.stalled[i] = 0
                self.LEF_birth(i)

    def _force_unload_and_rebirth(self, i: int):
        """Force-unload LEF i and immediately re-birth it."""
        self.occupied[self.LEFs[i, 0]] = 0
        self.occupied[self.LEFs[i, 1]] = 0
        self.stalled[i, :] = 0
        self.currentSpeeds[i] = 1.0
        self.phase[i] = 0
        self.dwell[i] = 0
        self.LEF_birth(i)

    def LEF_step(self):
        for i in range(self.numLEF):
            if self.mode == "asymmetric":
                if self.dwell[i] > self.dwell_time:
                    self._force_unload_and_rebirth(i)
                leg_id = self.directions[i]
                cur = self.LEFs[i, leg_id]
                stall_p = self.stallProbLeft[cur] if leg_id == 0 else self.stallProbRight[cur]
                
                # --- Core Hypothesis: Whether to enter phase II ---
                is_stable = (self.stalled[i, 0] == 1 and self.stalled[i, 1] == 1)
                
                if is_stable:
                    self.dwell[i] += 1
                
                # --- If colide with CTCF ---
                if np.random.random() < stall_p:
                    current_switch_p = self.switchProb * (self.stability_factor if is_stable else 1.0)
                    
                    if np.random.random() < current_switch_p:
                        self.directions[i] = 1 - leg_id 
                        leg_id = self.directions[i]
                        self.stalled[i, leg_id] = 0
                        self.currentSpeeds[i] = 1.0 if is_stable else self.speedMultiplier
                    else:
                        self.stalled[i, leg_id] = 1
                
                # --- Move ---
                if not self.stalled[i, leg_id]:
                    eff_speed = self.currentSpeeds[i]
                    if is_stable:
                        eff_speed = min(eff_speed, 1.0)
                    
                    attempts = int(eff_speed) + (1 if np.random.random() < (eff_speed % 1) else 0)
                    for _ in range(attempts):
                        cur = self.LEFs[i, leg_id]
                        if leg_id == 0:
                            if cur > 0 and self.occupied[cur - 1] == 0:
                                if np.random.random() > self.pause[cur]:
                                    self.occupied[cur-1], self.occupied[cur] = 1, 0
                                    self.LEFs[i, 0] = cur - 1
                        else:
                            if cur < self.numSite - 1 and self.occupied[cur + 1] == 0:
                                if np.random.random() > self.pause[cur]:
                                    self.occupied[cur+1], self.occupied[cur] = 1, 0
                                    self.LEFs[i, 1] = cur + 1

            if self.mode == "symmetric":
                for side in [0, 1]:
                    cur = self.LEFs[i, side]
                    stall_p = self.stallProbLeft[cur] if side == 0 else self.stallProbRight[cur]
                    if np.random.random() < stall_p:
                        self.stalled[i, side] = 1
                if np.random.random() < 0.5: # update for 1/2 speed
                    cur1, cur2 = self.LEFs[i]
                    if self.stalled[i, 0] == 0 and cur1 > 0 and self.occupied[cur1 - 1] == 0:
                        if np.random.random() > self.pause[cur1]:
                            self.occupied[cur1-1], self.occupied[cur1] = 1, 0
                            self.LEFs[i, 0] = cur1 - 1
                    if self.stalled[i, 1] == 0 and cur2 < self.numSite - 1 and self.occupied[cur2 + 1] == 0:
                        if np.random.random() > self.pause[cur2]:
                            self.occupied[cur2+1], self.occupied[cur2] = 1, 0
                            self.LEFs[i, 1] = cur2 + 1
    def step(self):
        self.LEF_death()
        self.LEF_step()

    def steps(self, N):
        for _ in range(N):
            self.step()

class LEFTranslocatorDynamicBoundary(LEFTranslocator):
    def __init__(self, numLEF, deathProb, stalledDeathProb, birthArray, pauseProb, 
                 stallProbLeft, stallProbRight, ctcfDeathProb, ctcfBirthProb, *args, 
                 mode="asymmetric", switchProb=0.0, speedMultiplier=1.0, stability_factor=0.1,
                 initalize_at_equilibrium_occupancy=True):
        self.ctcfDeathProb = ctcfDeathProb
        self.ctcfBirthProb = ctcfBirthProb
        self.stallProbLeft_init = np.copy(stallProbLeft)
        self.stallProbRight_init = np.copy(stallProbRight)
        
        super().__init__(numLEF, deathProb, stalledDeathProb, birthArray, pauseProb, 
                         stallProbLeft, stallProbRight, mode=mode, switchProb=switchProb, 
                         speedMultiplier=speedMultiplier, stability_factor=stability_factor)
        
        if initalize_at_equilibrium_occupancy:
            eq_occ = self.ctcfBirthProb / (self.ctcfBirthProb + self.ctcfDeathProb)
            self.stallProbLeft = (eq_occ > np.random.random(size=self.numSite)) * (self.stallProbLeft_init > 0)
            self.stallProbRight = (eq_occ > np.random.random(size=self.numSite)) * (self.stallProbRight_init > 0)

    def ctcf_death_left(self, i):
        self.stallProbLeft[i] = 0
        for j in range(self.numLEF):
            if i == self.LEFs[j, 0]: self.stalled[j, 0] = 0

    def ctcf_death_right(self, i):
        self.stallProbRight[i] = 0
        for j in range(self.numLEF):
            if i == self.LEFs[j, 1]: self.stalled[j, 1] = 0

    def ctcf_birth_left(self, i):
        self.stallProbLeft[i] = self.stallProbLeft_init[i]

    def ctcf_birth_right(self, i):
        self.stallProbRight[i] = self.stallProbRight_init[i]

    def step(self):
        updates = np.random.random((2, self.numSite))
        protected_idx = self.LEFs[self.dwell > 0].ravel()

        # left CTCF death
        mask_l = (updates[0] < self.ctcfDeathProb) & (self.stallProbLeft > 0)
        mask_l[protected_idx] = False 
        for ind in np.flatnonzero(mask_l):
            self.ctcf_death_left(ind)

        # right CTCF death
        mask_r = (updates[1] < self.ctcfDeathProb) & (self.stallProbRight > 0)
        mask_r[protected_idx] = False
        for ind in np.flatnonzero(mask_r):
            self.ctcf_death_right(ind)

        for ind in np.flatnonzero((updates[0] < self.ctcfBirthProb) * (self.stallProbLeft != self.stallProbLeft_init)):
            self.ctcf_birth_left(ind)
        for ind in np.flatnonzero((updates[1] < self.ctcfBirthProb) * (self.stallProbRight != self.stallProbRight_init)):
            self.ctcf_birth_right(ind)
        super().step()
