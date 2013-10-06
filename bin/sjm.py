class Job:

	time=None
	memory=None
	queue=None
	project=None
	status=None
	log_dir=None
	cmd_prefix=None
	cmd_separator='&&'
	name_prefix=None
	sge_options=None

	def __init__(self, name=None):
		self.name=name
		if self.name_prefix is not None and self.name is not None:
			self.name=self.name_prefix+self.name
		self.status=None
		self.cmds=[]
		self.dependents=[]

	def __str__(self):
		s='job_begin\n'
		if self.name is not None:
			s+='\tname %s\n'%self.name
		if self.time is not None:
			s+='\ttime %s\n'%self.time
		if self.memory is not None:
			s+='\tmemory %s\n'%self.memory
		if self.queue is not None:
			s+='\tqueue %s\n'%self.queue
		if self.project is not None:
			s+='\tproject %s\n'%self.project
		if self.status is not None:
			s+='\tstatus %s\n'%self.status
		if self.sge_options is not None:
			s+='\tsge_options %s\n'%self.sge_options
		if len(self.cmds)>0:
			s+='\tcmd_begin\n'
			s+=(' %s\n'%('' if self.cmd_separator is None else self.cmd_separator)).join(['\t\t%s %s'%(('' if self.cmd_prefix is None else self.cmd_prefix), cmd) for cmd in self.cmds])+"\n"
			s+='\tcmd_end\n'
		s+='job_end\n'
		return s

	def done(self):
		self.status='done'

	def append(self, cmd):
		self.cmds.append(cmd)
		return self

	def depend(self, *jobs):
		if jobs is not None:
			for job in jobs:
				if job is not None:
					self.dependents.append(job)
		return self

	def order(self, history=[]):
		s=''
		for dependent in self.dependents:
			s+=dependent.order(history)
			order=(dependent.name, self.name)
			if self.name is not None and order not in history:
				s+= "order %s before %s\n" % order
				history.append(order)
		return s

	def traverse(self, history=[]):
		s=''
                for dependent in self.dependents:
                        s+=dependent.traverse(history)
		if self.name is not None and self not in history:
			s+=str(self)
			history.append(self)
                return s

	def desc(self):
		s=self.traverse()
		s+=self.order()
		if self.log_dir is not None:
			s+='log_dir %s\n'%self.log_dir
		return s
