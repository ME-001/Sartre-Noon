<h3>
  Issues faced while working
</h3><l><hr><br>
<p>The Following are the issues that were faced while working on the project</p>
<ul>
  <li><b> 'STU' not defined:</b>b></li>
    <p>
      The "STU" parameter set that was used for running Sartre in upc mode needs to be manually added. 
    Tobias provided the required files both src and tables-related files, which needed to be placed in the right location.
    </p> 
    <p>
      The issue is that even after copying the files at the required location and then source building Sartre then also the issue persisted with <i>error Parameter set 'STU' not defined</i>.
    </p>
  <p>
    This was a Linux issue that if you download a file and copy it to another location then with that file an additional file generated at the new location named something 'zone change'.
    Remove these files and source build sartre.
    The Issue should be resolved.
  </p>
</ul>

**This is a markdown format file ".md" which doesn't need html formmating as far as I know.**
