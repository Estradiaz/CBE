

function grid(x_, y_, z_, id_){//(size_x_, size_y_, size_z_){
	//radiale darstellung von jedem punkt??
	//als knoten mit 8+6 nachbarelementen
	this.solid = true;
	
	//randevents laufen in den overflow
	//solution: fixed border with reflection -> needs no neighbours for any calculation
	// every grid element with any x,y,z == 0 is solid
	//ok now calculate correctly

	//this.drawSphere = function(r, d_){
		this.drawSphere = function(){
		var r = map(this.rho,1,3,0,2);
		var red = map(this.rho,0,255,0,2);
		push();
		translate(this.x*10, this.y*10, this.z*10);
		fill(red+80,0,100);
		sphere(r)//r);
		//console.log("r",r,"d",d_);
		pop();
	}

	
	this.getIDofNeighbour = function(dx_, dy_, dz_){
		var x = this.x + dx_; //integers!!!!
		var y = this.y + dy_;
		var z = this.z + dz_;

		return (x + (x - 1) * this.yn * this.zn + y*(1 + this.zn) - this.zn +  z);
	}

	this.setSizes = function(xn_, yn_, zn_){
		this.xn = xn_;
		this.yn = yn_;
		this.zn = zn_;
	}
	this.neighbours = [];
	this.pushNeighbours = function(grids){
		if(!this.solid)
		for (i = 1; i < this.e.length ; i++)  {
			if (this.e[i].x != 0 && this.e[i].y != 0 && this.e[i].z != 0){
				var idn = this.getIDofNeighbour(this.e[i].x,this.e[i].y,this.e[i].z)
				//console.log("inner grids",grids);
				//console.log("getIDofNeighbour", idn );
				//console.log("getsPushedas idn", grids[idn]);
				this.neighbours.push(grids[idn]);
				//console.log("this.neighbours", this.neighbours);
			}
		}
	}



	//relativbezug -> 
	this.x = x_;
	this.y = y_;
	this.z = z_;
	this.id = id_;

	this.dx = 1;
	
	this.dt = 1;
	
	this.c_max = 1; // CFL < 1 bedingung für stabilität euler verfaheren
	this.CFL = true;

	
	this.particles = []; 

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
	this.updatef = function(){
		this.f = this.f_new.slice(0);
	}
	this.f_new = [];
	this.w = [];
	this.algorithm_depth = 2; //set as default
	this.rho = random(2);//1;
		this.setRho = function(rho_){
			this.rho = rho_;
		}
		this.updateRho = function(){
			if (!this.n) this.set_n();
			var out = 0;
			for (var i = 0; i < this.n; i++) {
				out+= this.f[i];
			};
			this.rho = out;
			return this.rho;	
		}

	this.c2 = 1;
		grid.prototype.c2 = 1;
	this.u = createVector(0,0,0);
		this.setU = function(u_){
			this.u = u_.copy();
		}
		this.updateU = function(){
			if (!this.n) this.set_n();
			for (var i = 0; i < this.n; i++) {
				this.u.add(this.e[i].mult(this.f[i]));
			};
			this.CFL = (this.u.x+this.u.y+this.u.z < this.dx/this.dt);
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
	this.nnn = createVector(0,0,0);
	this.e.push(this.nnn);
	//Flächen
	this.pnn = createVector(1,0,0); //vorne
	this.e.push(this.pnn);
	this.npn = createVector(0,1,0); //links
	this.e.push(this.npn);
	this.nnp = createVector(0,0,1); //oben
	this.e.push(this.nnp);
	this.mnn = createVector(-1,0,0); //hinten
	this.e.push(this.mnn);
	this.nmn = createVector(0,-1,0); //rechts
	this.e.push(this.nmn);
	this.nnm = createVector(0,0,-1); //unten
	this.e.push(this.nnm);
	//Ecken
	this.pmm = createVector(1,-1,-1); //vorne_rechts_unten
	this.e.push(this.pmm);
	this.ppm = createVector(1,1,-1); //vorne_links_unten
	this.e.push(this.ppm);
	this.pmp= createVector(1,-1,1); //vorne_rechts_oben
	this.e.push(this.pmp);
	this.ppp= createVector(1,1,1); //vorne_links_oben
	this.e.push(this.ppp);
	this.mmm= createVector(-1,-1,-1); //hinten_rechts_unten
	this.e.push(this.mmm);
	this.mpm= createVector(-1,1,-1); //hinten_links_unten
	this.e.push(this.mpm);
	this.mmp= createVector(-1,-1,1); //hinten_rechts_oben
	this.e.push(this.mmp);
	this.mpp= createVector(-1,1,1); //hinten_links_oben
	this.e.push(this.mpp);
	this.set_n();
	}
	


	this.initialize_f = function(){
		//console.log(this.e);
		for (i = 0; i < this.e.length; i++){
			//grid.prototype.eqilibrium_distr = function(rho, wi, ei, u)
			var f_val = grid.prototype.eqilibrium_distr(this.rho, this.w[i], this.e[i], this.u, 1);
			this.f.push(f_val);
			this.f_new.push(f_val);
		}
	}
	// f_i(x+e*delta_t,t+delta_t) = 



	this.collision = function(){
		//update params
		this.updateRho();
		this.updateU();
		if(!this.n) this.set_n();
		for (var i = 0; i < this.n; i++) {
			this.f[i] = grid.prototype.bgk(this.f[i], this.tau, this.rho, this.w[i], this.e[i], this.u);
		};
	}

	//infinit space , think about boundaries !?!
	//O(n²)
	this.stream = function(){
		//check lfor solids and maybe use neighbours neighbour.x,y,z - this.x,y,z = ei(dx,dy,dz)
		// hashlist local id -> global id <- is in definition of neighbours[localId].id returns global
		// neighbours.getlocal(global) O(n) 
		// haslist gloabal id -> local id !!
//		for (var i = 0; i < n; i++) {
//			this.particles[i].f_new[grid.getOpositeIndexByDirection(this.e(i))] = max(0,this.f_i[i]);
//		};
		for ( i = 0 ; i < this.neighbours.length; i++){
			var vec = createVector(this.x - this.neighbours[i].x, this.y - this.neighbours[i].y, this.z - this.neighbours[i].z);
			var ivec = grid.prototype.getIndexByDirection(vec);
			if (!this.neighbours[i].solid){
				this.neighbours[i].f_new[grid.prototype.getOpositeIndexByDirection(vec)] = max(0,this.f_i[ivec]);
			}
			else{
				this.f_new(ivec) = -1 * max(0,this.f_i[ivec]) ;
			};
		};

	}
}

	grid.prototype.getIndexByDirection = function(Vector){
		var tmp
		
		for (var i = 0; i < this.n; i++) {
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

	grid.prototype.eqilibrium_distr = function(rho, wi, ei, u, c2_ ){
		//wi * rho * (1+ ( 3* dot(ei,u) / c2) + ...)
		if (rho === 0) return 0;
		var tmp = wi * rho;
		//console.log("e",ei);
		//console.log("u",u);
		var dots = p5.Vector.dot(ei,u);
		var sums1 = 1;
		var sums2 = dots;
		sums2 *= 3;
		
		sums2 /= c2_;
		var sums3 = dots*dots;
		sums3 /= c2_;
		sums3 /= c2_;
		sums3 *= 9;
		sums3 /= 2;
		var sums4 = 3/c2_;
		sums4 *= p5.Vector.dot(u,u);
		return tmp*(sums1+sums2+sums3+sums4);

	}

	grid.prototype.bgk = function(f, tau, rho, w, e, u){
		var	f_eq =grid.prototype.eqilibrium_distr(rho, w, e, u);
		var f_i = f;
		var F_ext_i = 0;
		return f_i + (1/tau) * (f_eq - f_i) - F_ext_i;
	}

grid.prototype.initializeGridArray = function(xn,yn,zn){
		var grids = [];
		var id = 0
		for (var i = 0; i <  xn; i++) {
			 for (var j = 0; j < yn; j++) {
			 	for (var k = 0; k < zn; k++) {
			 		var tmp = new grid(i,j,k,id);
			 		id += 1;
		 			tmp.solid = ( i == 0 || i == xn-1 || j == 0 || j == yn - 1 || k == 0 || k == zn - 1) ;
			 		tmp.setSizes(xn,yn,zn);
			 		tmp.inititalize_D3Q19();
			 		tmp.initialize_f();
			 		grids.push(tmp);
			 	};
			 };
		};
		for (var i = 0; i < id-1; i++) {
			//console.log("push",i,grids[i]);
			grids[i].pushNeighbours(grids);
		};
		return grids;
	};

grid.prototype.runInstance = function(grids){
	console.log("runInstance");
	var g_l = grids.length;
	for (var i = 0; i < g_l; i++){
		if (!grids[i].solid){
			grids[i].collision();
			grids[i].stream();

		}
	}
		for (var i = 0; i < g_l; i++){
		if (!grids[i].solid){
			grids[i].updatef();
			grids[i].updateRho();
			grids[i].updateU();		
		}
	}
	console.log("/runInstance");
}

//trommelwirbel .... 
