Main Methods
=================================

.. py:function:: pypka.getTitrableSites(pdb)

	Gets the all titrable sites from the pdb

   :param pdb: The filename a PDB file   
   :type pdb: string

   :return: All titrable sites residue numbers found in the pdb file
   :rtype: list


.. py:class:: pypka.Titration

	Main pypka class

	.. py:method:: getAverageProt(self, site, pH)
        
        Calculates the average protonation of a site at a given pH value

        :param site: Residue number of the site with the suffix 'N'
        			or 'C' to indicate a N-ter or C-ter, respectively.
        :type site: string

        :param pH: pH value
        :type pH: float 
        
        :returns: Average protonation of the site at the selected pH value
        :rtype: float

	.. py:method:: getProtStategetProtState(self, site, pH)
        
        Indicates the most probable protonation state of a site at a given pH value

       	:param site: Residue number of the site with the suffix 'N'
            or 'C' to indicate a N-ter or C-ter, respectively.
       	:type site:  str
        :param pH: pH value 
        :type pH:  float

        :returns: A tuple with two elements. Example: ('protonated', 0.9)

            The first element is a string indicating the most probable
            state of the site. It can be either 'protonated',
            'deprotonated' or 'undefined'. The 'undefined' state is
            prompted if the average protonation state is between 0.1 and 0.9.
          
            The second element is a float of the average protonation of the site
        :rtype: tuple

	.. py:method:: getParameters()

		Get the parameters used in the calculations

		:returns: Current state of DelPhi parameters
		:rtype: string