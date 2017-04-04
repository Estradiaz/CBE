function grid(x_, y_, z_){//(size_x_, size_y_, size_z_){
	//radiale darstellung von jedem punkt??
	//als knoten mit 8+6 nachbarelementen

	//relativbezug -> 
	this.x = x_;
	this.y = y_;
	this.z = z_;

	this.dx = 1;
	grid.prototype.dx = 1;
	this.dt = 1;
	grid.prototype.dt = 1;
	this.c_max = 1; // CFL < 1 bedingung für stabilität euler verfaheren
	this.CFL = true;

	//particle innerhalb des grid
	this.particles = []; //max particles = 1 -< fuse particles ~ higher density/pressure/mass
	//how to handle inner repulsive forces?
	//vielmehr strömung ~ impuls p = m * v(i[1to16])
	//1D:?
	this.globalgrid = [];
	this.n = 0;
		this.set_n = function (){
			this.n = this.e.length;
		}
		this.get_n = function(){
			return this.n;
		}
	this.e = [];
	this.f = [];
	this.f_new = [];
	this.w = [];
	this.algorithm_depth = 2; //set as default
	this.rho = 1;
		this.setRho = function(rho_){
			this.rho = rho_;
		}
		this.updateRho = function(){
			var out = 0;
			for (var i = 0; i < n; i++) {
				out+= this.f[i];
			};
			this.rho = out;
			return this.rho;	
		}

	this.c2 = 1;
		grid.prototype.c2 = 1;
	this.u = new p5.Vector(0,0,0);
		this.setU = function(u_){
			this.u = u_.copy();
		}
		this.updateU = function(){
			for (var i = 0; i < n; i++) {
				this.u.add(this.e[i].mult(f[i]));
			};
			this.CFL = (this.u.x+this.u.y+this.u.z < dx/dt);
			if (!this.CFL) console.log("CFL invalid: change dt- | dx+  ")
		}


	this.setDepth = function(new_depth_){
		algorithm_depth = new_depth_;
	}




	// 14 neighbours -> self:=(0,0,0)
	// (x,y,z) -> +|- permutation(1,0,0) -> +|- permutation(1,1,0) & (1,1,1)
	// notation: n := null, m:= minus 1, p:= plus 1
	this.inititalize_D3Q19 = function() {
	// as per D3Q15
	this.nnn = new p5.Vector(0,0,0);
	this.e.push(this.nnn);
	//Flächen
	this.pnn = new p5.Vector(1,0,0); //vorne
	this.e.push(this.pnn);
	this.npn = new p5.Vector(0,1,0); //links
	this.e.push(this.npn);
	this.nnp = new p5.Vector(0,0,1); //oben
	this.e.push(this.nnp);
	this.mnn = new p5.Vector(-1,0,0); //hinten
	this.e.push(this.mnn);
	this.nmn = new p5.Vector(0,-1,0); //rechts
	this.e.push(this.nmn);
	this.nnm = new p5.Vector(0,0,-1); //unten
	this.e.push(this.nnm);
	//Ecken
	this.pmm = new p5.Vector(1,-1,-1); //vorne_rechts_unten
	this.e.push(this.pmm);
	this.ppm = new p5.Vector(1,1,-1); //vorne_links_unten
	this.e.push(this.ppm);
	this.pmp= new p5.Vector(1,-1,1); //vorne_rechts_oben
	this.e.push(this.pmp);
	this.ppp= new p5.Vector(1,1,1); //vorne_links_oben
	this.e.push(this.ppp);
	this.mmm= new p5.Vector(-1,-1,-1); //hinten_rechts_unten
	this.e.push(this.mmm);
	this.mpm= new p5.Vector(-1,1,-1); //hinten_links_unten
	this.e.push(this.mpm);
	this.mmp= new p5.Vector(-1,-1,1); //hinten_rechts_oben
	this.e.push(this.mmp);
	this.mpp= new p5.Vector(-1,1,1); //hinten_links_oben
	this.e.push(this.mpp);
	this.set_n();
	}
	
	grid.prototype.getIndexByDirection = function(Vector){
		var tmp
		for (var i = 0; i < n; i++) {
			tmp = this.e[i].copy();
			if (tmp.normalize() == Vector.normalize()) return i;
		};

	}
		grid.prototype.getOpositeIndexByDirection = function(Vector){
		var tmp;
		for (var i = 0; i < n; i++) {
			tmp = this.e[i].copy();
			if (tmp.mult(-1).normalize() == Vector.normalize()) return i;
		};
		
	}

	grid.prototype.eqilibrium_distr = function(rho, wi, ei, u){
		//wi * rho * (1+ ( 3* dot(ei,u) / c2) + ...)
		if (rho === 0) return 0;
		var tmp = wi * rho;
		var dots = p5.Vector.dot(ei,u);
		var sums1 = 1;
		var sums2 = dots;
		sums2 *= 3;
		sums2 /= grid.prtotype.c2;
		var sums3 = dots*dots;
		sums3 /= grid.prtotype.c2;
		sums3 /= grid.prtotype.c2;
		sums3 *= 9;
		sums3 /= 2;
		var sums4 = 3/grid.prtotype.c2;
		sums4 *= p5.Vector.dot(u,u);
		return tmp*(sums1+sums2+sums3+sums4);

	}

	this.initialize_f = function(){
		var f_val = grid.eqilibrium_distr(this.rho, this.w[i], this.c2, this.e[i], this.u);
		this.f.push(f_val);
		this.f_new.push(f_val);
	}
	// f_i(x+e*delta_t,t+delta_t) = 
	grid.prototype.bgk = function(f, tau, rho, w, e, u){
		var	f_eq =grid.eqilibrium_distr(rho, w, e, u);
		var f_i = f;
		var F_ext_i = 0;
		return f_i + (1/tau) * (f_eq - f_i) - F_ext_i;
	}


	this.collision = function(){
		//update params
		this.updateRho();
		this.updateU();
		for (var i = 0; i < n; i++) {
			this.f[i] = grid.bgk(this.f[i], this.tau, this.rho, this.w[i], this.e[i], thiw.u);
		};	
	}

	//infinit space , think about boundaries !?!
	//O(n²)
	this.stream = function(){
		for (var i = 0; i < n; i++) {
			this.particles[i].f_new[grid.getOpositeIndexByDirection(this.e(i))] = max(0,this.f_i[i]);
		};
	}
}

